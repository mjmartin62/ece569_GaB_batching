###########################################################

## USER SPECIFIC DIRECTORIES ##

# CUDA directory:
#CUDA_ROOT_DIR=/usr/local/cuda

##########################################################

## CC COMPILER OPTIONS ##

# CC compiler options:
CC=g++
CC_FLAGS=
CC_LIBS=

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
#set(COMPUTE_CAPABILITY 60)
# set(CUDA_SEPARABLE_COMPILATION TRUE)
#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --gpu-architecture;sm_${COMPUTE_CAPABILITY}")
NVCC=nvcc
NVCC_FLAGS=-lm --gpu-architecture sm_60 -std=c++11 -Xcompiler -Wno-unused-function,-Wno-unused-parameter
NVCC_LIBS=

# CUDA library directory:
#CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
#CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= 

##########################################################

## Project file structure ##

# Source file directory:
USER    = $(whoami)
SRC_DIR = ~$(USER)/ece569/proj/source/GaB

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = ~$(USER)/ece569/proj/source/GaB

##########################################################

## Make variables ##

# Target executable name:
EXE = run_test

# Object files:
OBJS = $(OBJ_DIR)/main


##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
$(EXE) : $(OBJS)
	$(OBJS) $(SRC_DIR)/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 2 64 20 0
#	$(CC) $(CC_FLAGS) $(OBJS) -o $@ $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)

# Compile main .cpp file to object files:
#$(OBJ_DIR)/%.o : %.cpp
#	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile C++ source files to object files:
$(OBJ_DIR)/% : $(SRC_DIR)/main.cu 
	$(NVCC)  $<  $(SRC_DIR)/GaB.cu $(NVCC_FLAGS) -o $@ -I=$(INC_DIR) 


#	$(NVCC)  $<  $(SRC_DIR)/GaB.cu $(NVCC_FLAGS) -o run_test 
# Compile CUDA source files to object files:
#$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cu $(INC_DIR)/%.cuh
#$(OBJ_DIR)/%.o : $(SRC_DIR)/GaB.cu 
#	$(NVCC) $< $(NVCC_FLAGS) -o $@ $(NVCC_LIBS) -c

# Clean objects in object directory.
clean:
	$(RM) bin/* $(EXE)