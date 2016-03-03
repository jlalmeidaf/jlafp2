SDK_COMMON_DIR = /usr/local/cuda/include
SDK_COMMON_LIB_DIR = /usr/local/cuda
TOOLKIT_INC_DIR = /usr/local/cuda/include

ARCH = sm_20    #use sm_20 for Fermi, sm_12 for 1.2 and sm_13 for 1.3
CXX = nvcc
NVCC = nvcc

NVCCOPTIONS = -arch $(ARCH) #-use_fast_math 
NVCCFLAGS = -O3 $(NVCCOPTIONS) #-I $(TOOLKIT_INC_DIR) -I. -I $(SDK_COMMON_DIR)/inc/

CXXFLAGS = $(NVCCFLAGS)
NVCCLIBS = #-L $(SDK_COMMON_LIB_DIR)/lib -lcutil -lm

all: source directories

#------------------------------ seq.out --------------------------------
source: gpas

gpas: thread_manager.o \
	matches_manager.o \
	hi_res_timer.o \
	thread.o \
	needleman_wunsch_global.o \
	exceptions.o \
	mutex.o \
	main.o \
	sem.o \
	smith_waterman_local.o \
	arguments_manager.o \
	substitution_matrix.o \
	similarity_algorithm_cpu.o \
	needleman_wunsch_semiglobal.o \
	sequences.o \
	main_cu.o
		  $(NVCC) $(NVCCFLAGS) $(NVCCLIBS) thread_manager.o \
		matches_manager.o \
		hi_res_timer.o \
		thread.o \
		needleman_wunsch_global.o \
		exceptions.o \
		mutex.o \
		main.o \
		sem.o \
		smith_waterman_local.o \
		arguments_manager.o \
		substitution_matrix.o \
		similarity_algorithm_cpu.o \
		needleman_wunsch_semiglobal.o \
		sequences.o \
		main_cu.o\
		-o gpas


#------------------------------ CPP --------------------------------------------
thread_manager.o: thread_manager.cpp thread_manager.h thread.h \
 exceptions.h sem.h mutex.h
		$(CXX) $(CXXFLAGS)  -c ./thread_manager.cpp
matches_manager.o: matches_manager.cpp matches_manager.h exceptions.h
		$(CXX) $(CXXFLAGS)  -c ./matches_manager.cpp
hi_res_timer.o: hi_res_timer.cpp hi_res_timer.h
		$(CXX) $(CXXFLAGS)  -c ./hi_res_timer.cpp
thread.o: thread.cpp thread.h
		$(CXX) $(CXXFLAGS)  -c ./thread.cpp
needleman_wunsch_global.o: needleman_wunsch_global.cpp \
 needleman_wunsch_global.h back_up_struct.h exceptions.h sequences.h \
 substitution_matrix.h similarity_algorithm_cpu.h matches_manager.h \
 hi_res_timer.h
		$(CXX) $(CXXFLAGS)  -c ./needleman_wunsch_global.cpp
exceptions.o: exceptions.cpp exceptions.h
		$(CXX) $(CXXFLAGS)  -c ./exceptions.cpp
mutex.o: mutex.cpp mutex.h
		$(CXX) $(CXXFLAGS)  -c ./mutex.cpp
main.o: main.cpp sequences.h substitution_matrix.h exceptions.h \
 needleman_wunsch_global.h back_up_struct.h similarity_algorithm_cpu.h \
 matches_manager.h hi_res_timer.h needleman_wunsch_semiglobal.h \
 smith_waterman_local.h main_cu.h thread_manager.h thread.h sem.h mutex.h \
 arguments_manager.h
		$(CXX) $(CXXFLAGS)  -c ./main.cpp
sem.o: sem.cpp sem.h
		$(CXX) $(CXXFLAGS)  -c ./sem.cpp
smith_waterman_local.o: smith_waterman_local.cpp smith_waterman_local.h \
 back_up_struct.h exceptions.h sequences.h substitution_matrix.h \
 similarity_algorithm_cpu.h matches_manager.h hi_res_timer.h
		$(CXX) $(CXXFLAGS)  -c ./smith_waterman_local.cpp
arguments_manager.o: arguments_manager.cpp arguments_manager.h
		$(CXX) $(CXXFLAGS)  -c ./arguments_manager.cpp
substitution_matrix.o: substitution_matrix.cpp substitution_matrix.h \
 exceptions.h
		$(CXX) $(CXXFLAGS)  -c ./substitution_matrix.cpp
similarity_algorithm_cpu.o: similarity_algorithm_cpu.cpp \
 similarity_algorithm_cpu.h back_up_struct.h exceptions.h sequences.h \
 substitution_matrix.h matches_manager.h hi_res_timer.h
		$(CXX) $(CXXFLAGS)  -c ./similarity_algorithm_cpu.cpp
needleman_wunsch_semiglobal.o: needleman_wunsch_semiglobal.cpp \
 needleman_wunsch_semiglobal.h back_up_struct.h exceptions.h sequences.h \
 substitution_matrix.h similarity_algorithm_cpu.h matches_manager.h \
 hi_res_timer.h needleman_wunsch_global.h
		$(CXX) $(CXXFLAGS)  -c ./needleman_wunsch_semiglobal.cpp
sequences.o: sequences.cpp sequences.h substitution_matrix.h exceptions.h
		$(CXX) $(CXXFLAGS)  -c ./sequences.cpp
#------------------------------ CU ---------------------------------------------
main_cu.o : main_cu.cu \
    main_cu.h \
    sequences.h \
    substitution_matrix.h \
    exceptions.h \
    thread_manager.h \
    thread.h \
    sem.h \
    mutex.h \
    matches_manager.h \
    data_management.cu \
    data_management.cuh \
    cuda_declarations.h \
    needleman_wunsch_global_score_gpu.cu \
    alignment_score_gpu.cuh \
    hi_res_timer.h \
    needleman_wunsch_semiglobal_score_gpu.cu \
    smith_waterman_score_gpu.cu \
    needleman_wunsch_global_match_gpu.cu \
    alignment_match_gpu.cuh \
    needleman_wunsch_semiglobal_match_gpu.cu \
    smith_waterman_match_gpu.cu



	  $(NVCC) $(NVCCFLAGS) -c main_cu.cu -o main_cu.o
#------------------------------ DIRECTORIES ------------------------------------
directories: results 

results:
	 mkdir -p results
#-------------------------------- CLEANING -------------------------------------
clean: 
	rm -f thread_manager.o 
	rm -f matches_manager.o 
	rm -f hi_res_timer.o 
	rm -f thread.o 
	rm -f needleman_wunsch_global.o 
	rm -f exceptions.o 
	rm -f mutex.o 
	rm -f main.o 
	rm -f sem.o 
	rm -f smith_waterman_local.o 
	rm -f arguments_manager.o 
	rm -f substitution_matrix.o 
	rm -f similarity_algorithm_cpu.o 
	rm -f needleman_wunsch_semiglobal.o 
	rm -f sequences.o 
	rm -f gpas
	rm -f main_cu.o

#--------------------------------- OTHER ---------------------------------------
ptx:
	$(NVCC) $(NVCCFLAGS) -c main_cu.cu --ptxas-options=-v 

