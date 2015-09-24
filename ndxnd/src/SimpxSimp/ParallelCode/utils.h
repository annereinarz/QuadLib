void assert_correct_number_of_inputs (int actual, int expected) {
   if (actual-1 != expected){
      printf("Requires %d input arguments, got %d.\n", expected, actual-1);
      exit(1);
   }
}

int get_integer_argument(char* argumentString, char* errorMsg) {
   int argument;
   if (1 != sscanf(argumentString, "%d", &argument)) {
      printf("%s\n", errorMsg);
      exit(1);
   }
   return argument;
}
