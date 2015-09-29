"""
Created on Jul 19, 2013


This class represents the OptiType model for HLA typing based on NGS data

It is dependent on Coopr and uses an external ILP solver such as GLPK or CPLEX



@author: Benjamin Schubert
"""

from __future__ import division
from pyomo.environ import ConcreteModel, Set, Param, Var, Binary, Objective, Constraint, ConstraintList, maximize
from pyomo.opt import SolverFactory, TerminationCondition
from collections import defaultdict
import pandas as pd
import itertools


class OptiType(object):
    """
    classdocs

    """

    def __init__(self, cov, occ, groups_4digit, allele_table, beta, t_max_allele=2, solver="glpk", threads=1,
                 verbosity=0):
        """
        Constructor
        """

        self.__allele_table = allele_table
        self.__beta = float(beta)
        self.__t_max_allele = t_max_allele
        self.__solver = SolverFactory(solver)
        self.__threads = threads
        self.__opts = {"threads": threads} if threads > 1 else {}
        self.__verbosity = verbosity
        self.__changed = True  # model needs to know if it changed from last run or not
        self.__ks = 1
        self.__groups_4digit = groups_4digit

        loci_alleles = defaultdict(list)
        for type_4digit, group_alleles in groups_4digit.iteritems():
            # print type_4digit, group_alleles
            loci_alleles[type_4digit.split('*')[0]].extend(group_alleles)

        loci = loci_alleles

        self.__allele_to_4digit = {allele: type_4digit for type_4digit, group in groups_4digit.iteritems() for allele in
                                   group}

        '''
            generates the basic ILP model
        '''

        model = ConcreteModel()

        # init Sets
        model.LociNames = Set(initialize=loci.keys())
        model.Loci = Set(model.LociNames, initialize=lambda m, l: loci[l])

        L = list(itertools.chain(*loci.values()))
        reconst = {allele_id: 0.01 for allele_id in L if '_' in allele_id}
        R = set([r for (r, _) in cov.keys()])
        model.L = Set(initialize=L)
        model.R = Set(initialize=R)

        # init Params
        model.cov = Param(model.R, model.L, initialize=lambda model, r, a: cov.get((r, a), 0))
        model.reconst = Param(model.L, initialize=lambda model, a: reconst.get(a, 0))

        model.occ = Param(model.R, initialize=occ)
        model.t_allele = Param(initialize=self.__t_max_allele, mutable=True)

        model.beta = Param(initialize=self.__beta,
                           validate=lambda val, model: 0.0 <= float(self.__beta) <= 0.999,
                           mutable=True)
        model.nof_loci = Param(initialize=len(loci))

        # init variables
        model.x = Var(model.L, domain=Binary)
        model.y = Var(model.R, domain=Binary)

        model.re = Var(model.R, bounds=(0.0, None))
        model.hetero = Var(bounds=(0.0, model.nof_loci))

        # init objective
        model.read_cov = Objective(
            rule=lambda model: sum(model.occ[r] * (model.y[r] - model.beta * (model.re[r])) for r in model.R) - sum(
                model.reconst[a] * model.x[a] for a in model.L), sense=maximize)

        # init Constraints
        model.max_allel_selection = Constraint(model.LociNames, rule=lambda model, l: sum(
            model.x[a] for a in model.Loci[l]) <= model.t_allele)
        model.min_allel_selection = Constraint(model.LociNames,
                                               rule=lambda model, l: sum(model.x[a] for a in model.Loci[l]) >= 1)
        model.is_read_cov = Constraint(model.R,
                                       rule=lambda model, r: sum(model.cov[r, a] * model.x[a] for a in model.L) >=
                                                             model.y[r])
        model.heterozygot_count = Constraint(
            rule=lambda model: model.hetero >= sum(model.x[a] for a in model.L) - model.nof_loci)

        # regularization constraints
        model.reg1 = Constraint(model.R, rule=lambda model, r: model.re[r] <= model.nof_loci * model.y[r])
        model.reg2 = Constraint(model.R, rule=lambda model, r: model.re[r] <= model.hetero)
        model.reg3 = Constraint(model.R,
                                rule=lambda model, r: model.re[r] >= model.hetero - model.nof_loci * (1 - model.y[r]))

        # generate constraint list for solution enumeration
        model.c = ConstraintList()
        # Generate instance. Used to be .create() but deprecated since,
        # as ConcreteModels are instances on their own now.
        self.__instance = model

    def set_beta(self, beta):
        """
            Sets the parameter beta
        """
        self.__changed = True
        getattr(self.__instance, str(self.__instance.beta)).set_value(float(beta))

    def set_t_max_allele(self, t_max_allele):
        """
            Sets the upper bound of alleles selected per loci
        """
        self.__changed = True
        getattr(self.__instance, str(self.__instance.t_allele)).set_value(t_max_allele)

    def solve(self, ks):
        """
            solves the problem k times and discards the found solutions in the next run.
        """
        d = defaultdict(list)  # in there we store the typing +objective and generate afterwards a DatarFrame with it

        if self.__changed or self.__ks != ks:
            self.__ks = ks
            for k in xrange(ks):
                expr = 0

                self.__instance.preprocess()
                try:
                    res = self.__solver.solve(self.__instance, options=self.__opts, tee=self.__verbosity)
                except:
                    print ("WARNING: Solver does not support multi-threading. Please change the config"
                          " file accordingly. Falling back to single-threading.")
                    res = self.__solver.solve(self.__instance, options={}, tee=self.__verbosity)
                self.__instance.solutions.load_from(res)  # solution loading changed recently.

                # if self.__verbosity > 0:
                #     res.write(num=1)

                if res.solver.termination_condition != TerminationCondition.optimal:
                    print "Optimal solution hasn't been obtained. This is a terminal problem."  # TODO message, exit
                    break

                selected = []
                indices = []
                encountered_4digit = []
                for j in self.__instance.x:
                    if self.__allele_to_4digit[j][0] in 'HJG':
                        if 0.99 <= self.__instance.x[j].value <= 1.01:
                            selected.append(j)
                        indices.append(j)
                        continue
                    if 0.99 <= self.__instance.x[j].value <= 1.01:
                        selected.append(j)
                        exp_i = 0
                        exp_i += self.__instance.x[j]
                        if self.__allele_to_4digit[j] in encountered_4digit:
                            continue
                        encountered_4digit.append(self.__allele_to_4digit[j])
                        for i_allele in self.__groups_4digit[self.__allele_to_4digit[j]]:
                            if self.__instance.x[i_allele].value <= 0:
                                exp_i += self.__instance.x[i_allele]
                            indices.append(i_allele)
                        expr += (1.0 - exp_i)
                zero_indices = set([j for j in self.__instance.x]).difference(set(indices))
                for j in zero_indices:
                    expr += self.__instance.x[j]

                self.__instance.c.add(expr >= 1)

                # if self.__verbosity > 0:
                #     print selected
                #     self.__instance.c.pprint()
                aas = [self.__allele_to_4digit[x].split('*')[0] for x in selected]
                c = dict.fromkeys(aas, 1)
                for i in xrange(len(aas)):
                    if aas.count(aas[i]) < 2:
                        d[aas[i] + "1"].append(selected[i])
                        d[aas[i] + "2"].append(selected[i])
                    else:
                        d[aas[i] + str(c[aas[i]])].append(selected[i])
                        c[aas[i]] += 1

                nof_reads = sum((self.__instance.occ[j] * self.__instance.y[j].value for j in self.__instance.y))
                # if self.__verbosity > 0:
                #     print "Obj", res.Solution.Objective.__default_objective__.Value
                d['obj'].append(self.__instance.read_cov())
                d['nof_reads'].append(nof_reads)

            self.__instance.c.clear()
            self.__changed = False
            self.__enumeration = pd.DataFrame(d)

            # self.__rank()
            return self.__enumeration
        else:
            return self.__enumeration

    def solve_for_k_alleles(self, k, ks=1):
        """
            EXPERIMENTAL!

            generates a solution without the regularization term and only k selected alleles
        """
        if k < int(self.__instance.nof_loci.value) or k > int(self.__instance.nof_loci.value * self.__t_max_allele):
            raise Warning("k " + str(k) + " is out of range [" + str(self.__instance.nof_loci.value) + "," + str(
                self.__instance.nof_loci * self.__t_max_allele) + "]")


        # copy the instance
        inst = self.__instance.clone()
        # set beta = 0  # because we do homozygosity calling manually
        getattr(inst, str(inst.beta)).set_value(float(0.0))

        inst.del_component("heterozygot_count")
        inst.del_component("reg1")
        inst.del_component("reg2")
        inst.del_component("reg3")
        # generate constraint which allows only k alleles to be selected
        expr1 = 0
        for j in inst.x:
            expr1 += inst.x[j]

        inst.c.add(expr1 == k)
        d = defaultdict(list)

        for _ in xrange(ks):
            inst.preprocess()
            try:
                res = self.__solver.solve(inst, options=self.__opts, tee=self.__verbosity)
            except:
                print ("WARNING: Solver does not support multi-threading. Please change the config"
                      " file accordingly. Falling back to single-threading.")
                res = self.__solver.solve(inst, options={}, tee=self.__verbosity)
            inst.solutions.load_from(res)

            if self.__verbosity > 0:
                res.write(num=1)

            if res.solver.termination_condition != TerminationCondition.optimal:
                print "Optimal solution hasn't been obtained. This is a terminal problem."  # TODO message, exit
                break

            selected = []
            expr = 0

            indices = []
            encountered_4digit = []
            for j in inst.x:
                if 0.99 <= inst.x[j].value <= 1.01:
                    exp_i = 0
                    selected.append(j)
                    exp_i += inst.x[j]
                    if self.__allele_to_4digit[j] in encountered_4digit:
                        continue

                    encountered_4digit.append(self.__allele_to_4digit[j])
                    for i_allele in self.__groups_4digit[self.__allele_to_4digit[j]]:
                        if inst.x[i_allele].value <= 0:
                            exp_i += inst.x[i_allele]
                        indices.append(i_allele)
                    expr += (1 - exp_i)
            zero_indices = set([j for j in inst.x]).difference(set(indices))
            for j in zero_indices:
                expr += inst.x[j]

            inst.c.add(expr >= 1)

            if self.__verbosity > 0:
                print selected
            aas = [self.__allele_to_4digit[x].split('*')[0] for x in selected]
            c = dict.fromkeys(aas, 1)
            for i in xrange(len(aas)):
                if aas.count(aas[i]) < 2:
                    d[aas[i] + "1"].append(selected[i])
                    d[aas[i] + "2"].append(selected[i])
                else:
                    d[aas[i] + str(c[aas[i]])].append(selected[i])
                    c[aas[i]] += 1

            # print "Obj", res.Solution.Objective.__default_objective__.Value
            nof_reads = sum((inst.occ[j] * inst.y[j].value for j in inst.y))
            d['obj'].append(inst.read_cov())
            d['nof_reads'].append(nof_reads)

        return pd.DataFrame(d)

    def solve_fixed_typing(self, fixed_alleles):
        """
            EXPERIMENTAL!

            forces the allele to pic a 4-digit of the provided alleles
        """
        k = len(set(fixed_alleles))
        if k < int(self.__instance.nof_loci.value) or k > int(self.__instance.nof_loci.value * self.__t_max_allele):
            raise Warning("k " + str(k) + " is out of range [" + str(self.__instance.nof_loci.value) + "," + str(
                self.__instance.nof_loci * self.__t_max_allele) + "]")


        # copy the instance
        inst = self.__instance.clone()
        # set beta = 0 because we do homozygocity calling manually
        getattr(inst, str(inst.beta)).set_value(float(0.0))

        inst.del_component("heterozygot_count")
        inst.del_component("reg1")
        inst.del_component("reg2")
        inst.del_component("reg3")
        # generate constraint which allows only k alleles to be selected
        expr1 = 0
        for j in inst.x:
            expr1 += inst.x[j]
        inst.c.add(expr1 == k)

        # generate for each of the provided alleles the fixation constraint:
        for a in set(fixed_alleles):
            expr_f = 0
            print self.__groups_4digit
            for ids in self.__groups_4digit[a]:
                print ids
                expr_f += inst.x[ids]
            inst.c.add(expr_f == 1)

        d = defaultdict(list)

        inst.preprocess()
        try:
            res = self.__solver.solve(inst, options=self.__opts, tee=self.__verbosity)
        except:
            print ("WARNING: Solver does not support multi-threading. Please change the config"
                  " file accordingly. Falling back to single-threading.")
            res = self.__solver.solve(inst, options={}, tee=self.__verbosity)
        inst.solutions.load_from(res)

        opt_ids = [j for j in inst.x if 0.99 <= inst.x[j].value <= 1.01]

        aas = [self.__allele_to_4digit[x].split('*')[0] for x in opt_ids]
        c = dict.fromkeys(aas, 1)
        for i in xrange(len(aas)):
            if aas.count(aas[i]) < 2:
                d[aas[i] + "1"].append(opt_ids[i])
                d[aas[i] + "2"].append(opt_ids[i])
            else:
                d[aas[i] + str(c[aas[i]])].append(opt_ids[i])
                c[aas[i]] += 1
        nof_reads = sum((inst.occ[j] * inst.y[j].value for j in inst.y))
        d['obj'].append(self.inst.read_cov())
        d['nof_reads'].append(nof_reads)

        return pd.DataFrame(d)

    def enumerate_allele_wise(self):
        """
            EXPERIMENTAL!

            fixes all but one allele and solves it again to investigate the influence of this
            particular allele on the objective value.
        """
        d = defaultdict(list)

        self.__instance.preprocess()
        try:
            res = self.__solver.solve(self.__instance, options=self.__opts, tee=self.__verbosity)
        except:
            print ("WARNING: Solver does not support multi-threading. Please change the config"
                  " file accordingly. Falling back to single-threading.")
            res = self.__solver.solve(self.__instance, options={}, tee=self.__verbosity)
        self.__instance.solutions.load_from(res)

        opt_ids = [j for j in self.__instance.x if 0.99 <= self.__instance.x[j].value <= 1.01]

        aas = [self.__allele_to_4digit[x].split('*')[0] for x in opt_ids]
        c = dict.fromkeys(aas, 1)
        for i in xrange(len(aas)):
            if aas.count(aas[i]) < 2:
                d[aas[i] + "1"].append(opt_ids[i])
                d[aas[i] + "2"].append(opt_ids[i])
            else:
                d[aas[i] + str(c[aas[i]])].append(opt_ids[i])
                c[aas[i]] += 1
        nof_reads = sum((self.__instance.occ[j] * self.__instance.y[j].value for j in self.__instance.y))
        d['obj'].append(self.__instance.read_cov())
        d['nof_reads'].append(nof_reads)
        d['discarded'].append(0)

        for j in opt_ids:
            if self.__verbosity > 0:
                self.__instance.c.pprint()
            self.__instance.c.clear()
            # fix all but j'th variable
            fix = 0
            for i in opt_ids:
                if i != j:
                    fix += (1 - self.__instance.x[i])
            self.__instance.c.add(fix == 0.0)

            # discard j'th allele and all its 4digit equivalent alleles form the next solution
            discard = 0
            for k in self.__groups_4digit[self.__allele_to_4digit[j]]:
                discard += self.__instance.x[k]
            self.__instance.c.add(discard == 0.0)

            # solve with new constraints
            self.__instance.preprocess()
            try:
                res = self.__solver.solve(self.__instance, tee=self.__verbosity)  # ,tee=True) verbose solvinf
                self.__instance.solutions.load_from(res)
            except:
                print Warning("There is no replacement for allele " + self.__allele_to_4digit[j])
                continue

            selected = [al for al in self.__instance.x if 0.99 <= self.__instance.x[al].value <= 1.01]
            aas = [self.__allele_to_4digit[x].split('*')[0] for x in selected]
            c = dict.fromkeys(aas, 1)
            for q in xrange(len(aas)):
                if aas.count(aas[q]) < 2:
                    d[aas[q] + "1"].append(selected[q])
                    d[aas[q] + "2"].append(selected[q])
                else:
                    d[aas[q] + str(c[aas[q]])].append(selected[q])
                    c[aas[q]] += 1
            nof_reads = sum((self.__instance.occ[h] * self.__instance.y[h].value for h in self.__instance.y))
            d['obj'].append(self.__instance.read_cov())
            d['nof_reads'].append(nof_reads)
            d['discarded'].append(j)
        return pd.DataFrame(d)

    def solve_enforced_zygosity(self, gosity_dict):
        """
            EXPERIMENTAL!

            solves the ilp without regularization but enforced homo/heterozygosity for each locus
            @param gosity_dict: a dictionary with all loci as keys and value = number of alleles per locus (default is 2)
        """

        inst = self.__instance.clone()
        # set beta = 0 because we do homozygocity calling manually
        getattr(inst, str(inst.beta)).set_value(float(0.0))

        inst.del_component("heterozygot_count")
        inst.del_component("reg1")
        inst.del_component("reg2")
        inst.del_component("reg3")

        # now delete max_allele_constraint and reconstruct it again
        inst.del_component("max_allel_selection")
        for locus in inst.LociNames:
            cons = 0
            for a in inst.Loci[locus]:
                cons += inst.x[a]
            inst.c.add(cons <= gosity_dict.get(locus, 2))

        d = defaultdict(list)

        inst.preprocess()
        try:
            res = self.__solver.solve(inst, options=self.__opts, tee=self.__verbosity)
        except:
            print ("WARNING: Solver does not support multi-threading. Please change the config"
                  " file accordingly. Falling back to single-threading.")
            res = self.__solver.solve(inst, options={}, tee=self.__verbosity)
        inst.solutions.load_from(res)

        selected = [al for al in inst.x if 0.99 <= inst.x[al].value <= 1.01]
        aas = [self.__allele_to_4digit[x].split('*')[0] for x in selected]
        c = dict.fromkeys(aas, 1)
        for q in xrange(len(aas)):
            if aas.count(aas[q]) < 2:
                d[aas[q] + "1"].append(selected[q])
                d[aas[q] + "2"].append(selected[q])
            else:
                d[aas[q] + str(c[aas[q]])].append(selected[q])
                c[aas[q]] += 1
        nof_reads = sum((inst.occ[h] * inst.y[h].value for h in inst.y))
        d['obj'].append(inst.read_cov())
        d['nof_reads'].append(nof_reads)
        return pd.DataFrame(d)
