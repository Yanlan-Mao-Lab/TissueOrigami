################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../SourceCode/Node.cpp \
../SourceCode/Prism.cpp \
../SourceCode/ReferenceShapeBase.cpp \
../SourceCode/ShapeBase.cpp \
../SourceCode/Simulation.cpp \
../SourceCode/main.cpp 

OBJS += \
./SourceCode/Node.o \
./SourceCode/Prism.o \
./SourceCode/ReferenceShapeBase.o \
./SourceCode/ShapeBase.o \
./SourceCode/Simulation.o \
./SourceCode/main.o 

CPP_DEPS += \
./SourceCode/Node.d \
./SourceCode/Prism.d \
./SourceCode/ReferenceShapeBase.d \
./SourceCode/ShapeBase.d \
./SourceCode/Simulation.d \
./SourceCode/main.d 


# Each subdirectory must supply rules for building sources it contributes
SourceCode/%.o: ../SourceCode/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


