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
CND_CONF=Debug_score
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/bic_scoring_function.o \
	${OBJECTDIR}/top_p_constraint.o \
	${OBJECTDIR}/score_calculator.o \
	${OBJECTDIR}/ad_tree.o \
	${OBJECTDIR}/bayesian_network.o \
	${OBJECTDIR}/bdeu_scoring_function.o \
	${OBJECTDIR}/ad_node.o \
	${OBJECTDIR}/fnml_scoring_function.o \
	${OBJECTDIR}/sparse_parent_tree.o \
	${OBJECTDIR}/vary_node.o \
	${OBJECTDIR}/file_pattern_database.o \
	${OBJECTDIR}/score_cache.o \
	${OBJECTDIR}/log_likelihood_calculator.o \
	${OBJECTDIR}/tarjans_algorithm.o \
	${OBJECTDIR}/ad_tree_scoring_main.o \
	${OBJECTDIR}/dynamic_pattern_database.o \
	${OBJECTDIR}/static_pattern_database.o \
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
LDLIBSOPTIONS=-L${HOME}/local/lib -lboost_system -lboost_program_options -lboost_timer -lrt -lboost_chrono -lboost_thread

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/score-debug

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/score-debug: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/score-debug ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/bic_scoring_function.o: bic_scoring_function.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/bic_scoring_function.o bic_scoring_function.cpp

${OBJECTDIR}/dynamic_pattern_database.o: dynamic_pattern_database.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -s -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/dynamic_pattern_database.o dynamic_pattern_database.cpp

${OBJECTDIR}/static_pattern_database.o: static_pattern_database.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -Wall -s -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/static_pattern_database.o static_pattern_database.cpp

${OBJECTDIR}/top_p_constraint.o: top_p_constraint.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/top_p_constraint.o top_p_constraint.cpp

${OBJECTDIR}/score_calculator.o: score_calculator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/score_calculator.o score_calculator.cpp

${OBJECTDIR}/ad_tree.o: ad_tree.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ad_tree.o ad_tree.cpp

${OBJECTDIR}/bayesian_network.o: bayesian_network.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/bayesian_network.o bayesian_network.cpp

${OBJECTDIR}/bdeu_scoring_function.o: bdeu_scoring_function.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/bdeu_scoring_function.o bdeu_scoring_function.cpp

${OBJECTDIR}/ad_node.o: ad_node.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ad_node.o ad_node.cpp

${OBJECTDIR}/fnml_scoring_function.o: fnml_scoring_function.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/fnml_scoring_function.o fnml_scoring_function.cpp

${OBJECTDIR}/sparse_parent_tree.o: sparse_parent_tree.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/sparse_parent_tree.o sparse_parent_tree.cpp

${OBJECTDIR}/vary_node.o: vary_node.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/vary_node.o vary_node.cpp

${OBJECTDIR}/file_pattern_database.o: file_pattern_database.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/file_pattern_database.o file_pattern_database.cpp

${OBJECTDIR}/score_cache.o: score_cache.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/score_cache.o score_cache.cpp

${OBJECTDIR}/log_likelihood_calculator.o: log_likelihood_calculator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/log_likelihood_calculator.o log_likelihood_calculator.cpp

${OBJECTDIR}/tarjans_algorithm.o: tarjans_algorithm.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/tarjans_algorithm.o tarjans_algorithm.cpp

${OBJECTDIR}/ad_tree_scoring_main.o: ad_tree_scoring_main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ad_tree_scoring_main.o ad_tree_scoring_main.cpp

${OBJECTDIR}/hugin_structure_writer.o: hugin_structure_writer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I${HOME}/local/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/hugin_structure_writer.o hugin_structure_writer.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/score-debug

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
