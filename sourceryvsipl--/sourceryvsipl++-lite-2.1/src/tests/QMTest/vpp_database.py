########################################################################
#
# File:   vpp_database.py
# Author: Stefan Seefeld
# Date:   2005-04-15
#
# Contents:
#   VPPTest, VPPDatabase
#
# Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. 
#
########################################################################

########################################################################
# Imports
########################################################################

import fnmatch
import os
import qm
import qm.test.base
from   qm.fields import *
from   qm.executable import *
from   qm.test.resource import Resource
from   qm.test.database import *
from   qm.test.classes.compilation_test import *
import dircache

########################################################################
# Classes
########################################################################

class ParallelService(Resource):

    def SetUp(self, context, result):

        setup = Executable()
        command = []
        self.halt_command = []
        
        service = context.get('par_service.name')
        command = context.get('par_service.boot', '').split()
        self.halt_command = context.get('par_service.halt', '').split()

        if command:
            status = setup.Run(command)
            result.CheckExitStatus('ParallelService', ' '.join(command), status)

    def CleanUp(self, result):

        if self.halt_command:
            command = self.halt_command
            cleanup = Executable()
            status = cleanup.Run(command)
            result.CheckExitStatus('ParallelService', ' '.join(command), status)
        
            
class VPPDatabase(Database):
    """A 'VPPDatabase' stores the vsipl++ regression tests."""

    srcdir = TextField(title = "Source Directory",
                       description = "The root of the test suite's source tree.")
    _is_generic_database = True
    

    def __init__(self, path, arguments):

        self.label_class = "file_label.FileLabel"
        self.modifiable = "false"
        # Initialize the base class.
        super(VPPDatabase, self).__init__(path, arguments)

        
    def GetRoot(self):

        return self.srcdir


    def GetSubdirectories(self, directory):

        dirname = os.path.join(self.GetRoot(), directory)
        return [subdir for subdir in dircache.listdir(dirname)
                if (os.path.isdir(os.path.join(dirname, subdir)) and
                    subdir not in ('data', 'CVS', '.svn', 'QMTest'))]


    def GetIds(self, kind, directory = "", scan_subdirs = 1):

        dirname = os.path.join(self.GetRoot(), directory)
        if not os.path.isdir(dirname):
            raise NoSuchSuiteError, directory

        if kind == Database.TEST:
            datadir = os.path.join(dirname, 'data')
            if os.path.isdir(datadir):
                ids = [self.JoinLabels(directory, f)
                       for f in dircache.listdir(datadir)
                       if os.path.isfile(os.path.join(datadir, f))]
            else:
                ids = [self.JoinLabels(directory, f)
                       for f in dircache.listdir(dirname)
                       if os.path.splitext(f)[1] == '.cpp'] 
        elif kind == Database.RESOURCE:
            datadir = os.path.join(dirname, 'data')
            if os.path.isdir(datadir):
                ids = [self.JoinLabels(directory, f)
                       for f in dircache.listdir(dirname)
                       if os.path.splitext(f)[1] == '.cpp'] 
            else:
                ids = []
            
        else: # SUITE
            ids = [self.JoinLabels(directory, d)
                   for d in self.GetSubdirectories(directory)]

        if scan_subdirs:
            for subdir in dircache.listdir(dirname):
                if (subdir != 'data'
                    and os.path.isdir(os.path.join(dirname, subdir))):
                    dir = self.JoinLabels(directory, subdir)
                    ids.extend(self.GetIds(kind, dir, True))
        return ids
    

    def GetExtension(self, item_id):

        if not item_id:
            return DirectorySuite(self, item_id)
            
        if item_id == 'parallel_service':
            return ParallelService({},
                                   qmtest_id = item_id, qmtest_database = self)

        id_components = self.GetLabelComponents(item_id)
        # 'data' subdirectories have special meaning, and so
        # are not allowed as label components.
        if 'data' in id_components:
            return None

        resources = ['parallel_service']

        dirname = os.path.join(self.GetRoot(), *id_components[:-1])
        basename = id_components[-1]

        if os.path.isdir(os.path.join(dirname, basename)):
            return DirectorySuite(self, item_id)

        # If <dirname>/data is an existing directory...
        if os.path.isdir(os.path.join(dirname, 'data')):
            # ...and <basename> ends with '.cpp',...
            if os.path.splitext(basename)[1] == '.cpp':
                
                # ...<dirname>/<basename> is a resource.
                path = os.path.join(dirname, basename)
                arguments = {}
                arguments['source_files'] = [path]
                arguments['executable'] = os.path.splitext(basename)[0]
                arguments['resources'] = resources
                return CompiledResource(arguments,
                                        qmtest_id = item_id,
                                        qmtest_database = self)
            else:
                # ...<dirname>/<basename> is a test.
                path = os.path.join(dirname, 'data', basename)
                if not os.path.isfile(path):
                    return None
            
                src = [f for f in dircache.listdir(dirname)
                       if os.path.splitext(f)[1] == '.cpp']
                # There must be exactly one source file, which
                # is our resource.
                if len(src) > 1:
                    raise DatabaseError,\
                          'multiple source files found in %s'%dirname
                resources.append(self.JoinLabels(*(id_components[:-1] + src)))

                arguments = {}
                arguments['args'] = [path]
                arguments['resources'] = resources

                return ExecutableTest(arguments,
                                      qmtest_id = item_id,
                                      qmtest_database = self)
        
        else:
            path = os.path.join(dirname, basename)
            if os.path.isfile(path) and os.path.splitext(path)[1] == '.cpp':

                # <path> is a test.
            
                arguments = {}
                arguments['source_files'] = [path]
                arguments['executable'] = os.path.splitext(basename)[0]
                arguments['resources'] = resources

                return CompilationTest(arguments,
                                       qmtest_id = item_id,
                                       qmtest_database = self)

            elif os.path.isdir(path):
                return DirectorySuite(self, item_id)
            
        return None
