# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from math import sin, cos, tan, log, log10, e, pi
from operator import add, sub, mul, truediv, pow
from pyparsing import Literal, CaselessLiteral, Word, Combine, Optional, ZeroOrMore, Forward, nums, alphas


class Eval(object):
    def __init__(self, expression):
        self.__expr_stack = self.__parser(expression)

    def _init_unpickle(self, config):
        self.__expr_stack = config

    def pickle(self):
        return self.__expr_stack.copy()

    @classmethod
    def unpickle(cls, config):
        if not isinstance(config, list):
            raise Exception('Invalid config')
        obj = cls.__new__(cls)
        obj._init_unpickle(config.copy())
        return obj

    def calc(self, val):
        return self.__evaluate_stack([str(val) if x == 'X' else x for x in self.__expr_stack])

    @classmethod
    def __parser(cls, expression):
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

