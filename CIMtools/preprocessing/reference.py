# -*- coding: utf-8 -*-
#
#  Copyright 2016-2018 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from math import sin, cos, tan, log, log10, e, pi
from operator import add, sub, mul, truediv, pow
from pandas import DataFrame, Index, read_csv
from pyparsing import Literal, CaselessLiteral, Word, Combine, Optional, ZeroOrMore, Forward, nums, alphas
from CGRtools.containers import ReactionContainer, MoleculeContainer, CGRContainer
from sklearn.base import BaseEstimator
from .common import TransformerMixin
from ..exceptions import ConfigurationError


class MetaReference(BaseEstimator, TransformerMixin):
    def __init__(self, data):
        self.data = data
        self.__init()

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_MetaReference__')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__init()

    def set_params(self, **params):
        if params:
            super().set_params(**params)
            self.__init()
        return self

    def __init(self):
        try:
            self.__ext_header = self.__prepare_ext_header(self.data)
            self.__extension = self.__prepare_ext(self.data)
        except Exception as ex:
            raise ConfigurationError(ex)

    def transform(self, x, return_domain=False):
        x = super().transform(x)

        res = []
        for s in x:
            tmp = {}
            res.append(tmp)
            for key, value in s.meta.items():
                if key in self.__extension:
                    ek = self.__extension[key]
                    if callable(ek):
                        tmp[key] = ek(value)
                    else:
                        tmp.update(ek[value])

        res = DataFrame(res, columns=self.__ext_header, index=Index(range(len(res)), name='structure'))
        if return_domain:
            return x, ~res.isnull().any(axis=1)
        return res

    @staticmethod
    def __prepare_ext_header(data):
        """
        :param data: dict
        :return: list of strings. descriptors header
        """
        tmp = []
        for orig_key in sorted(data):
            operation = data[orig_key]
            if isinstance(operation, dict):
                tmp.extend(sorted(list(operation.values())[0]))  # get columns for replacement
            else:
                tmp.append(orig_key)
        return tmp

    @staticmethod
    def __prepare_ext(data):
        tmp = {}
        for k, v in data.items():
            if isinstance(v, dict):
                tmp[k] = v
            elif v:
                tmp[k] = Eval(v)
            else:
                tmp[k] = float
        return tmp

    _dtype = ReactionContainer, MoleculeContainer, CGRContainer


class Eval:
    def __init__(self, expression):
        self.__expr_stack = self.__parser(expression)

    def __call__(self, value):
        return self.__evaluate_stack([str(value) if x == 'X' else x for x in self.__expr_stack])

    @staticmethod
    def __parser(expression):
        """ adopted from Paul McGuire example. http://pyparsing.wikispaces.com/file/view/fourFn.py
        """
        expr_stack = []

        def push_first(strg, loc, toks):
            expr_stack.append(toks[0])

        def push_u_minus(strg, loc, toks):
            if toks and toks[0] == '-':
                expr_stack.append('unary -')

        point = Literal('.')
        _e = CaselessLiteral('E')
        fnumber = Combine(Word('+-' + nums, nums) +
                          Optional(point + Optional(Word(nums))) +
                          Optional(_e + Word('+-' + nums, nums)))
        ident = Word(alphas, alphas + nums + '_$')

        plus = Literal("+")
        minus = Literal("-")
        mult = Literal("*")
        div = Literal("/")
        lpar = Literal("(").suppress()
        rpar = Literal(")").suppress()
        addop = plus | minus
        multop = mult | div
        expop = Literal("^")
        _pi = CaselessLiteral("PI")
        x = CaselessLiteral("X")

        expr = Forward()
        atom = (Optional("-") + (x | _pi | _e | fnumber | ident + lpar + expr + rpar).setParseAction(push_first) |
                (lpar + expr.suppress() + rpar)).setParseAction(push_u_minus)

        factor = Forward()
        factor << atom + ZeroOrMore((expop + factor).setParseAction(push_first))

        term = factor + ZeroOrMore((multop + factor).setParseAction(push_first))
        expr << term + ZeroOrMore((addop + term).setParseAction(push_first))

        expr.parseString(expression)
        return expr_stack

    @classmethod
    def __evaluate_stack(cls, s):
        op = s.pop()
        if op == 'unary -':
            return -cls.__evaluate_stack(s)
        if op in '+-*/^':
            op2 = cls.__evaluate_stack(s)
            op1 = cls.__evaluate_stack(s)
            return cls.__opn[op](op1, op2)
        elif op == "PI":
            return pi
        elif op == "E":
            return e
        elif op in cls.__fn:
            return cls.__fn[op](cls.__evaluate_stack(s))
        elif op[0].isalpha():
            return 0
        else:
            return float(op)

    __opn = {'+': add, '-': sub, '*': mul, '/': truediv, '^': pow}

    __fn = dict(sin=sin, cos=cos, tan=tan, lg=log10, ln=log, abs=abs, trunc=lambda a: int(a), round=round,
                sgn=lambda a: (1 if a > 0 else -1) if abs(a) > 1e-12 else 0)


def prepare_metareference(params, csv='EXTKEY'):
    """
    :param params: list of key[:value] strings
    key is metadata key for calculation.
    if value skipped then metadata key will be converted to float
    if value started with '=' then metadata will be evaluated by equation presented after =
    else value is path to CSV file with header. one of the header's column should be named EXTKEY
    or any other setted by csv kwarg.
    this column wll be used as keys for metadata replacement by values from other columns
    :param csv: header column used as keys

    :return: MetaReference instance
    """
    extdata = {}
    for p in params:
        if ':' in p:
            ext, val = p.split(':', maxsplit=1)
            if val:
                if val[0] == '=':
                    extdata[ext] = val[1:]
                else:
                    tmp = read_csv(val)
                    k = tmp.pop(csv)
                    v = tmp.rename(columns=lambda x: '%s.%s' % (ext, x))
                    extdata[ext] = {x: y.to_dict() for x, (_, y) in zip(k, v.iterrows())}
            else:
                extdata[ext] = None
        else:
            extdata[p] = None
    return MetaReference(extdata)


__all__ = ['MetaReference', 'prepare_metareference']
