################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../SourceCode/Analysis.cpp \
../SourceCode/CellMigration.cpp \
../SourceCode/GrowthFunctionBase.cpp \
../SourceCode/Lumen.cpp \
../SourceCode/ModelInputObject.cpp \
../SourceCode/NewtonRaphsonSolver.cpp \
../SourceCode/Node.cpp \
../SourceCode/Prism.cpp \
../SourceCode/RandomGenerator.cpp \
../SourceCode/ReferenceShapeBase.cpp \
../SourceCode/ShapeBase.cpp \
../SourceCode/Simulation.cpp \
../SourceCode/main.cpp 

OBJS += \
./SourceCode/Analysis.o \
./SourceCode/CellMigration.o \
./SourceCode/GrowthFunctionBase.o \
./SourceCode/Lumen.o \
./SourceCode/ModelInputObject.o \
./SourceCode/NewtonRaphsonSolver.o \
./SourceCode/Node.o \
./SourceCode/Prism.o \
./SourceCode/RandomGenerator.o \
./SourceCode/ReferenceShapeBase.o \
./SourceCode/ShapeBase.o \
./SourceCode/Simulation.o \
./SourceCode/main.o 

CPP_DEPS += \
./SourceCode/Analysis.d \
./SourceCode/CellMigration.d \
./SourceCode/GrowthFunctionBase.d \
./SourceCode/Lumen.d \
./SourceCode/ModelInputObject.d \
./SourceCode/NewtonRaphsonSolver.d \
./SourceCode/Node.d \
./SourceCode/Prism.d \
./SourceCode/RandomGenerator.d \
./SourceCode/ReferenceShapeBase.d \
./SourceCode/ShapeBase.d \
./SourceCode/Simulation.d \
./SourceCode/main.d 


# Each subdirectory must supply rules for building sources it contributes
SourceCode/%.o: ../SourceCode/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	gccg++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


