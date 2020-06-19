#
#    GPT - Grid Python Toolkit
#    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)
#                        Mattia Bruno
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
from urllib import request
import gpt, os

class repo:
    def __init__(self,first,second=None):
        self.files = {}
        self.baseurl = first
        self.git_branch = second
        self.cache_dir = '/tmp/__gptrepositorycache__'
        
    def addfile(self,fname):
        fn=fname.split('/')[-1]
        self.files[fn]=fname
        
    def ls(self):
        ls = ', '.join([k for k in self.files.keys()])
        gpt.message(f'Repository files: {ls}')
    
    def cached(self):
        if not os.path.exists(self.cache_dir):
            gpt.message('Nothing cached')
            return
        n=0
        for fn in self.files:
            fname=f'{self.cache_dir}/{fn}'
            if (gpt.rank()==0) and (os.path.exists(fname)):
                gpt.message(f'File {fn} cached')
                n+=1
        if (gpt.rank()==0) and (n==0):
            gpt.message('Nothing cached')
            
    def clean_cache(self):
        if not os.path.exists(self.cache_dir):
            return
        for f in os.listdir(self.cache_dir):
            os.remove(f)
        os.rmdir(self.cached_dir)
        gpt.message('Repository cache cleaned')
        
    def download(self,dest=None):
        if dest is None:
            dest = self.cache_dir
        os.makedirs(self.cache_dir,exist_ok=True)
        for fn in self.files:
            if self.git_branch is None:
                url=f'{self.baseurl}/{self.files[fn]}'
            else:
                url=f'{self.baseurl}/raw/{self.git_branch}/{self.files[fn]}'                
            fname=f'{dest}/{fn}'
            if (gpt.rank()==0) and (not os.path.exists(fname)):
                dt=-gpt.time()
                fname, header = request.urlretrieve(url, filename=fname)
                dt+=gpt.time()
                MB=os.path.getsize(fname) / 1024.**2
                gpt.message('Downloaded file %s at %g MB/s' % (fn,MB/dt))
        gpt.barrier()
        
    def remove(self,dest=None):
        if dest is None:
            dest = self.cache_dir
        for fn in self.files:
            dst=fn.split('/')[-1]
            fname=f'{dest}/{dst}'
            if (gpt.rank()==0) and (os.path.exists(fname)):
                os.remove(fname)
                gpt.message('Removed file %s' % dst)
        gpt.barrier()
