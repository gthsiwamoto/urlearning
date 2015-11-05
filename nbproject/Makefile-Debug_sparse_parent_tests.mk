#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug_sparse_parent_tests
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/top_p_constraint.o \
	${OBJECTDIR}/priority_queue.o \
	${OBJECTDIR}/markov_network_score_calculator.o \
	${OBJECTDIR}/sparse_parent_bitwise.o \
	${OBJECTDIR}/bayesian_network.o \
	${OBJECTDIR}/sparse_parent_tests_main.o \
	${OBJECTDIR}/static_pattern_database.o \
	${OBJECTDIR}/dynamic_pattern_database.o \
	${OBJECTDIR}/ad_tree.o \
	${OBJECTDIR}/sparse_parent_tree.o \
	${OBJECTDIR}/markov_network_bdeu_scoring_function.o \
	${OBJECTDIR}/file_pattern_database.o \
	${OBJECTDIR}/score_cache.o \
	${OBJECTDIR}/tarjans_algorithm.o \
	${OBJECTDIR}/sparse_parent_list.o \
	${OBJECTDIR}/hugin_structure_writer.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-std=c++0x
CXXFLAGS=-std=c++0x

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L${HOME}/local/lib -lboost_timer -lboost_program_options -lboost_system -lboost_chrono -lboost_thread

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/sparse_parent_tests-debug

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/sparse_parent_tests-debug: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/sparse_parent_tests-debug ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/top_p_constraint.o: top_p_constraint.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/top_p_constraint.o top_p_constraint.cpp

${OBJECTDIR}/priority_queue.o: priority_queue.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/priority_queue.o priority_queue.cpp

${OBJECTDIR}/markov_network_score_calculator.o: markov_network_score_calculator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/markov_network_score_calculator.o markov_network_score_calculator.cpp

${OBJECTDIR}/sparse_parent_bitwise.o: sparse_parent_bitwise.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/sparse_parent_bitwise.o sparse_parent_bitwise.cpp

${OBJECTDIR}/bayesian_network.o: bayesian_network.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/bayesian_network.o bayesian_network.cpp

${OBJECTDIR}/sparse_parent_tests_main.o: sparse_parent_tests_main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/sparse_parent_tests_main.o sparse_parent_tests_main.cpp

${OBJECTDIR}/static_pattern_database.o: static_pattern_database.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/static_pattern_database.o static_pattern_database.cpp

${OBJECTDIR}/dynamic_pattern_database.o: dynamic_pattern_database.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/dynamic_pattern_database.o dynamic_pattern_database.cpp

${OBJECTDIR}/ad_tree.o: ad_tree.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ad_tree.o ad_tree.cpp

${OBJECTDIR}/sparse_parent_tree.o: sparse_parent_tree.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/sparse_parent_tree.o sparse_parent_tree.cpp

${OBJECTDIR}/markov_network_bdeu_scoring_function.o: markov_network_bdeu_scoring_function.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/markov_network_bdeu_scoring_function.o markov_network_bdeu_scoring_function.cpp

${OBJECTDIR}/file_pattern_database.o: file_pattern_database.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/file_pattern_database.o file_pattern_database.cpp

${OBJECTDIR}/score_cache.o: score_cache.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/score_cache.o score_cache.cpp

${OBJECTDIR}/tarjans_algorithm.o: tarjans_algorithm.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/tarjans_algorithm.o tarjans_algorithm.cpp

${OBJECTDIR}/sparse_parent_list.o: sparse_parent_list.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/sparse_parent_list.o sparse_parent_list.cpp

${OBJECTDIR}/hugin_structure_writer.o: hugin_structure_writer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/hugin_structure_writer.o hugin_structure_writer.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/sparse_parent_tests-debug

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
