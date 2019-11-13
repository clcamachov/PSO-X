################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/config.cpp \
../src/hcjob.cpp \
../src/particle.cpp \
../src/problem.cpp \
../src/pso.cpp \
../src/rng.cpp \
../src/swarm.cpp \
../src/utils.cpp 

OBJS += \
./src/config.o \
./src/hcjob.o \
./src/particle.o \
./src/problem.o \
./src/pso.o \
./src/rng.o \
./src/swarm.o \
./src/utils.o 

CPP_DEPS += \
./src/config.d \
./src/hcjob.d \
./src/particle.d \
./src/problem.d \
./src/pso.d \
./src/rng.d \
./src/swarm.d \
./src/utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/gsl/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


