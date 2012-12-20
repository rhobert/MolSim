################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/test/LinkedCellParticleContainerTest.cpp \
../src/test/ParticleContainerTest.cpp \
../src/test/SimpleParticleContainerTest.cpp \
../src/test/TestSettings.cpp 

OBJS += \
./src/test/LinkedCellParticleContainerTest.o \
./src/test/ParticleContainerTest.o \
./src/test/SimpleParticleContainerTest.o \
./src/test/TestSettings.o 

CPP_DEPS += \
./src/test/LinkedCellParticleContainerTest.d \
./src/test/ParticleContainerTest.d \
./src/test/SimpleParticleContainerTest.d \
./src/test/TestSettings.d 


# Each subdirectory must supply rules for building sources it contributes
src/test/%.o: ../src/test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/chris/workspace/MolSim/libxsd -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


