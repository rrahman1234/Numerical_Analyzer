
CC=g++


#EXTRA_LIB_FLAGS = -larmadillo
#LIB_FLAGS = $(EXTRA_LIB_FLAGS)
INCLIB=numerical_solver.h newton_raphson_solve.h quadratic_solve.h bisection_solve.h linear_solve.h sparse_solve.h secant_solve.h 

# "If any of second_inheritance.o pet.o animal.o human.o have changed
# since second_i changed, then rebuild second_i ... and the following
# lines tell you how to do it!"

main_i:       main.o numerical_solver.o newton_raphson_solve.o quadratic_solve.o bisection_solve.o linear_solve.o sparse_solve.o secant_solve.o
		$(CC) -o main_i main.o numerical_solver.o newton_raphson_solve.o quadratic_solve.o bisection_solve.o linear_solve.o sparse_solve.o secant_solve.o -larmadillo 

# Standard rules ... because we have four compiles all of which can be
# done with the same template / pattern

%.o:            %.cpp $(INCLIB)
		$(CC) -std=c++14 -I/Users/rrahman/Documents/Research_Learn_Stuffs/C++_Train/coding_exercise/coding_exercise/2021_drill_exercises/lib_tests/eigen_libTest/eigen/  -c -o $@ $<

# This "make test" is a bit informal - but it does show people how to
# run the code and tells them the results to expect.
# $(CC) -std=c++14 -I/Users/rrahman/Documents/Research_Learn_Stuffs/C++_Train/2021_drill_exercises/lib_tests/eigen_libTest/eigen/ --verbose -c -o $@ $<

test:
	@./main_i
	@echo "Test Completed - should have give 4 line report"

clean:
	@rm -f *.o
	@echo "All cleaned!"
