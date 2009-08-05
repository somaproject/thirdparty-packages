########################################################################
#
# File:   database.py
# Author: Stefan Seefeld
# Date:   2008-05-03
#
# Contents:
#   Test, Resource, and Database classes for the Sourcery VSIPL++ test suite.
#
# Copyright (c) 2008 by CodeSourcery, Inc.  All rights reserved. 
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
from   qm.test.classes.compilation_test import ExecutableTest
from   qm.test.classes.compilation_test_database import *
from   qm.test.directory_suite import DirectorySuite
import dircache

########################################################################
# Classes
########################################################################

class CompiledResource(Resource):
    """A CompiledResource fetches compilation parameters from environment
    variables CPPFLAGS, <lang>_options, and <lang>_ldflags in addition
    to the CompilerTable-related parameters."""

    options = SetField(TextField(), computed="true")
    ldflags = SetField(TextField(), computed="true")
    source_files = SetField(TextField(), computed="true")
    executable = TextField(computed="true")
    language = TextField()

    def SetUp(self, context, result):

        self._context = context
        self._compiler = CompilationTest({'options':self.options,
                                          'ldflags':self.ldflags,
                                          'source_files':self.source_files,
                                          'executable':self.executable,
                                          'language':self.language,
                                          'execute':False},
                                         qmtest_id = self.GetId(),
                                         qmtest_database = self.GetDatabase())
        
        self._compiler.Run(context, result)
        directory = self._compiler._GetDirectory(context)
        self._executable = os.path.join(directory, self.executable)
        context['CompiledResource.executable'] = self._executable
        

    def CleanUp(self, result):

        self._compiler._RemoveDirectory(self._context, result)


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
        
            
class Database(CompilationTestDatabase):
    """'Database' stores the Sourcery VSIPL++ test suite.

    In addition to the CompilationTestDatabase behavior, we must:

    * make all tests depend on the ParallelService resource
    * add special handling for directories containing 'data/' subdirs.


    """
    enable_cvsip_tests = TextField()
    excluded_subdirs = SetField(TextField(),
                                default_value = ['QMTest', '.svn', 'data', 'build'],
                                description="Subdirectories not to scan for tests.",
                                computed="true")


    def GetSubdirectories(self, directory):

        subdirs = super(Database, self).GetSubdirectories(directory)
        if not directory and self.enable_cvsip_tests != '1':
            if 'cvsip' in subdirs: del subdirs[subdirs.index('cvsip')]
        return subdirs


    def GetIds(self, kind, directory = '', scan_subdirs = 1):

        dirname = os.path.join(self.srcdir, directory)
        if os.path.isdir(os.path.join(dirname, 'data')):
            if kind == Database.TEST:
                return [self.JoinLabels(directory, f)
                        for f in dircache.listdir(os.path.join(dirname, 'data'))
                        if f not in self.excluded_subdirs]
            else:
                return []
        else:
            return super(Database, self).GetIds(kind, directory, scan_subdirs)

        if not os.path.isdir(dirname):
            raise NoSuchSuiteError, directory

        if kind == Database.TEST:
            ids = [self.JoinLabels(directory, f)
                   for f in dircache.listdir(dirname)
                   if (os.path.isfile(os.path.join(dirname, f)) and
                       os.path.splitext(f)[1] in self.test_extensions)]

        elif kind == Database.RESOURCE:
            ids = []
            
        else: # SUITE
            ids = [self.JoinLabels(directory, d)
                   for d in self.GetSubdirectories(directory)
                   if d not in self.excluded_subdirs]

        if scan_subdirs:
            for subdir in dircache.listdir(dirname):
                if (subdir not in self.excluded_subdirs
                    and os.path.isdir(os.path.join(dirname, subdir))):
                    dir = self.JoinLabels(directory, subdir)
                    ids.extend(self.GetIds(kind, dir, True))

        return ids
    

    def GetExtension(self, id):

        if not id:
            return DirectorySuite(self, id)
            
        elif id == 'compiler_table':
            return CompilerTable({}, qmtest_id = id, qmtest_database = self)

        elif id == 'parallel_service':
            return ParallelService({}, qmtest_id = id, qmtest_database = self)

        resources = ['compiler_table', 'parallel_service']

        id_components = self.GetLabelComponents(id)
        # 'data' subdirectories have special meaning, and so
        # are not allowed as label components.
        if 'data' in id_components:
            return None

        dirname = os.path.join(self.srcdir, *id_components[:-1])
        basename = id_components[-1]

        file_ext = os.path.splitext(basename)[1]

        # If <dirname>/data is an existing directory...
        if os.path.isdir(os.path.join(dirname, 'data')):

            if file_ext in self.test_extensions:

                executable = os.path.splitext(os.path.basename(id))[0]
                if sys.platform == 'win32':
                    executable += '.exe'

                # ...<dirname>/<basename> is a resource.
                src = os.path.abspath(os.path.join(self.srcdir, id))
                return self._MakeItem(id,
                                      CompiledResource,
                                      language=self.test_extensions[file_ext],
                                      source_files=[src],
                                      executable=executable,
                                      resources=resources)
            else:
                # ...<dirname>/<basename> is a test.
                path = os.path.join(dirname, 'data', basename)
                if not os.path.isfile(path):
                    return None

                src = [f for f in dircache.listdir(dirname)
                       if os.path.splitext(f)[1] in self.test_extensions]
                # There must be exactly one source file, which
                # is our resource.
                if len(src) > 1:
                    raise DatabaseError('multiple source files found in %s'%dirname)

                resources.append(self.JoinLabels(*(id_components[:-1] + src)))
                return self._MakeItem(id,
                                      ExecutableTest,
                                      resources=resources,
                                      args=[path])
            

        src = os.path.join(self.srcdir, id)
        if file_ext in self.test_extensions and os.path.isfile(src):
            executable = os.path.splitext(os.path.basename(id))[0]
            if sys.platform == 'win32':
                executable += '.exe'

            return self._MakeItem(id,
                                  CompilationTest,
                                  language=self.test_extensions[file_ext],
                                  source_files=[src],
                                  executable=executable,
                                  resources=resources)

        elif os.path.isdir(src):
            if not basename in self.excluded_subdirs:
                return DirectorySuite(self, id)

        else:
            return None


    def _MakeItem(self, id, class_, **args):

        return class_(args, qmtest_id = id, qmtest_database = self)
