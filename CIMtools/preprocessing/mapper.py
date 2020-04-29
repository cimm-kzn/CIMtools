# -*- coding: utf-8 -*-
#
#  Copyright 2015-2019 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools.containers import ReactionContainer
from CGRtools import RDFRead, RDFWrite
from os import walk, path
from pathlib import Path
from shutil import rmtree
from sklearn.base import BaseEstimator, TransformerMixin
from subprocess import run
from tempfile import mkdtemp
from ..exceptions import ConfigurationError
from ..utils import iter2array


class Mapper(BaseEstimator,TransformerMixin):
    def __init__(self):
        pass
    
    def transform(self,x):
        x = iter2array(x,dtype=(ReactionContainer))
        
        work_dir = Path(mkdtemp(prefix='mapp_'),dir=str(Path.home())) 
        inp_file = work_dir / 're_map.rdf'
        out_folder = work_dir / 'results'
    
        for dirpath, dirnames, filenames in walk("."):
            for filename in filenames:
                if filename == 'aam-utils-j8 (1).jar':
                    path_to_jar = path.join(dirpath, filename)
        
        with RDFWrite(inp_file) as f:
            for num,r in enumerate(x):
                r.meta['Id'] = num
                f.write(r)
        try:
            p = run(['java','-jar',path_to_jar,'-j','MAPPING','-i',inp_file,'-o',out_folder,'-rdf_id','Id','-min','-mixture'])    
        except FileNotFoundError as e: 
            raise ConfigurationError from e

        if p.returncode != 0: 
            raise ConfigurationError('execution failed') 
        
        x_out = RDFRead(Path(out_folder/'MAX_reactions.rdf')).read()
        if len(x) != len(x_out):
            raise ValueError('invalid data')
        
        rmtree(work_dir)
        return x_out
               
__all__ = ['Mapper']    
