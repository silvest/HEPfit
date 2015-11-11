# check python version
import sys
if sys.version_info.major  < 3:
    print('Please use python3!')
    exit(1)

# command-line argument
if (len(sys.argv) != 3):
    print ('Usage: python3', sys.argv[0], '<project name> <configurations.xml>')
    exit(1)
project = sys.argv[1]
confFile = sys.argv[2]

#####################################################################

# parse configurations.xml
import xml.etree.ElementTree as ET
tree = ET.parse(confFile)
root = tree.getroot()

# fetch release configuration
confs = root.find('confs')
for conf in confs:
    if conf.get('name') == 'Release':
        ReleaseConf = conf

# fetch include-paths 
IncPaths = []
incDir = ReleaseConf.find('compileType/ccTool/incDir')
if incDir != None:
    for pElem in incDir.iter('pElem'):
        if 'BAT' in pElem.text:
            IncPaths.append('${BAT_INCLUDE_PATH}')
        elif 'FeynHiggs' in pElem.text:
            IncPaths.append('${FH_INC}')
        elif 'LoopTools' in pElem.text:
            IncPaths.append('${LOOPTOOLS_DIR}/include')
        #elif 'boost' in pElem.text:
        #    IncPaths.append('${BOOST_INC}')
        else:
            if project != 'Analysis':
                IncPaths.append('${CMAKE_CURRENT_SOURCE_DIR}/' + pElem.text)
IncPaths.append('${BOOST_INC}')
                
# fetch command-line options for cctool
ccTool_commandLine = ReleaseConf.find('compileType/ccTool/commandLine') 
if ccTool_commandLine != None:
    cxxOpts = ccTool_commandLine.text

# fetch lib paths and libraries, linker options for Analysis 
if project == 'Analysis':
    LibPaths = []
    Libs = [] 
    LibModules = [] 
    linkerAddLib = ReleaseConf.find('compileType/linkerTool/linkerAddLib')
    linkerLibItems = ReleaseConf.find('compileType/linkerTool/linkerLibItems')
    linker_commandLine = ReleaseConf.find('compileType/linkerTool/commandLine')
    if linker_commandLine != None:
        linkerOpts = linker_commandLine.text
    for pElem in linkerAddLib.iter('pElem'):
        if 'BAT' in pElem.text:
            LibPaths.append('${BAT_LIB}')
        elif 'FeynHiggs/lib64' in pElem.text:
            LibPaths.append('${FH_LIB}')
        elif 'LoopTools' in pElem.text:
            LibPaths.append('${LT_LIB}')
        elif 'boost' in pElem.text:
            LibPaths.append('${BOOST_LIB}')
        else:
            if 'LoopTools' not in pElem.text and 'FeynHiggs/lib' not in pElem.text:
                LibPaths.append('${CMAKE_CURRENT_SOURCE_DIR}/' + pElem.text)
    for libProject in linkerLibItems.iter('linkerLibProjectItem'):
        #LibPaths.append('${CMAKE_CURRENT_SOURCE_DIR}/' + libProject.find('makeArtifact').get('PL'))
        ##
        #lib = libProject.find('makeArtifact').get('OP')
        #lib = lib.replace('${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/', '')
        ##
        lib = libProject.find('makeArtifact').get('PL')
        lib = lib.replace('../', '')
        LibModules.append(lib)
    for lib in linkerLibItems.iter('linkerLibLibItem'):
        Libs.append(lib.text)

# fetch header and source files
Headers = []
Srcs = []
for item in ReleaseConf.iter('item'):
    path = item.get('path')
    tool = item.get('tool')
    if 'src' in path or 'main.cpp' in path: # exclude test files
        if '3' in tool:
            Headers.append(path)
        elif '1' in tool:
            Srcs.append(path)

#####################################################################

print('cmake_minimum_required(VERSION 2.6)\n', sep='')

# project name
print('project(', project, ')\n', sep='')

# set Release as default
print('if(NOT CMAKE_BUILD_TYPE)')
print('  set(CMAKE_BUILD_TYPE Release)')
print('endif(NOT CMAKE_BUILD_TYPE)\n')

# cxx flags
#if ccTool_commandLine != None:
#    print('set(CMAKE_CXX_FLAGS \"-Wall\")', sep='') 
#    cxxOptList = cxxOpts.split('`')
#    for cxxOpt in cxxOptList:
#        if 'gsl-config' in cxxOpt:
#            print('set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} ${GSL_CFLAGS}\")', sep='') 
#        elif 'root-config' in cxxOpt:
#            print('set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS}\")', sep='') 
#        elif cxxOpt != '' and cxxOpt != ' ':
#            print('set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}', cxxOpt, '\")', sep='') 
print('set(CMAKE_CXX_FLAGS \"-Wall ${GSL_CFLAGS} ${ROOT_CFLAGS} -DBOOST_NO_CXX11_RVALUE_REFERENCES\")', sep='') 
print()

# include paths
if IncPaths:
    if project == 'Analysis':
        print('foreach(_project ${PROJECTLIST} MonteCarlo)', sep='')
        print('  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../${_project}/src)', sep='')
        print('endforeach(_project)', sep='')
    for incPath in IncPaths:
        if incPath == '${LT_INC}':
            print('if(LOOPTOOLS)', sep='')
            print('  include_directories(${LT_INC})', sep='')
            print('endif(LOOPTOOLS)', sep='')
        else:
            print('include_directories(', incPath, ')', sep='')
    print()

# linker flags for Analysis
if project == 'Analysis':
    if linker_commandLine != None:
        print('set(CMAKE_CXX_LINK_FLAGS \"\")', sep='') 
        linkerOptList = linkerOpts.split('`')
        for linkerOpt in linkerOptList:
            if 'gsl-config --libs' in linkerOpt:
                #print('set(CMAKE_CXX_LINK_FLAGS \"${CMAKE_CXX_LINK_FLAGS} ${GSL_LIBS}\")', sep='') 
                Libs.append('${GSL_LIBS}')
            elif 'root-config --libs' in linkerOpt:
                #print('set(CMAKE_CXX_LINK_FLAGS \"${CMAKE_CXX_LINK_FLAGS} ${ROOT_LIBS}\")', sep='') 
                Libs.append('${ROOT_LIBS}')
            elif 'root-config --ldflags' in linkerOpt:
                #print('set(CMAKE_CXX_LINK_FLAGS \"${CMAKE_CXX_LINK_FLAGS} ${ROOT_LDFLAGS}\")', sep='') 
                Libs.append('${ROOT_LDFLAGS}')
            elif '-lMinuit' in linkerOpt:
                Libs.append('-lMinuit')
            elif linkerOpt != '' and linkerOpt != ' ':
                print('set(CMAKE_CXX_LINK_FLAGS \"${CMAKE_CXX_LINK_FLAGS}', linkerOpt, '\")', sep='')
        print()

# library paths for Analysis
if project == 'Analysis':
    if LibPaths != None:
        for libPath in LibPaths:
            if libPath == '${LT_LIB}':
                print('if(LOOPTOOLS)')
                print('  link_directories(${LT_LIB})')
                print('endif(LOOPTOOLS)')
            else:
                print('link_directories(', libPath, ')', sep='')
        print()

# target executable/library and source files
if project == 'Analysis':
    print('if(${LIBTYPE} STREQUAL \"OBJECT\")')
    print('  foreach(_project MonteCarlo ${PROJECTLIST})')
    print('    set(target_obj_list ${target_obj_list} $<TARGET_OBJECTS:${_project}>)')
    print('  endforeach(_project)')
    print('  add_executable(analysis', end='')
    for src in Srcs:
        print(' ', src, sep='', end='')
    #for lib in LibModules:
    #    print(' $<TARGET_OBJECTS:', lib, '>', sep='', end='') # requires cmake > v2.8.8
    print(' ${target_obj_list})')

    print('else()')

    print('  add_executable(analysis', end='')
    for src in Srcs:
        print(' ', src, sep='', end='')
    print(')')
    print('  target_link_libraries(analysis', end='')
    print(' MonteCarlo ${PROJECTLIST}', end='')
    #for lib in LibModules:
    #    print(' ', lib, sep='', end='')
    print(')')
    print('endif()')
    print()

    for lib in Libs:
        if lib == 'gfortran':
            print('target_link_libraries(analysis ${LFORTRAN})')

        ## use absolute paths to the libraries
        #elif lib == 'BATmodels':
        #    print('target_link_libraries(analysis ${LBATM})')
        #elif lib ==  'BAT':
        #    print('target_link_libraries(analysis ${LBAT})')
        #elif lib == 'FH':
        #    print('target_link_libraries(analysis ${LFH})')
        #elif lib == 'ooptools':
        #    print('target_link_libraries(analysis ${LLT})')
        #elif lib == 'boost_program_options':
        #    print('target_link_libraries(analysis ${LBOOST})')

        ## use -looptools
        elif lib == 'ooptools':
            print('target_link_libraries(analysis ${LIBOOPTOOLS})')

        else:
            print('target_link_libraries(analysis ', lib, ')', sep='')
    print()
    print('add_dependencies(analysis BATBUILD)')
    print('add_dependencies(analysis FHBUILD)')
    print('if(MPIBAT)')
    print('  target_link_libraries(analysis ${MPI_LIBRARIES})')
    print('endif()')
    print('if(MPI_COMPILE_FLAGS)')
    print('  set_target_properties(analysis PROPERTIES COMPILE_FLAGS ${MPI_COMPILE_FLAGS})')
    print('endif()')
    print('if(MPI_LINK_FLAGS)')
    print('  set_target_properties(analysis PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})')
    print('endif()')
    print()
    print('INSTALL(TARGETS analysis DESTINATION bin COMPONENT executable)')

else:
    print('file(GLOB srcs \"${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp\")')
    print('add_library(', project, ' ${LIBTYPE}', sep='', end='')
    print(' ${srcs}', sep='', end='')
    #for src in Srcs:
    #    print(' ', src, sep='', end='')
    print(')')
    print()
    print('file(GLOB headers \"${CMAKE_CURRENT_SOURCE_DIR}/src/*.h\")')
    print('INSTALL(FILES ${headers} DESTINATION include/HEPfit COMPONENT header)')
print()



