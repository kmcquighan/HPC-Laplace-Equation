void  itoa    (int i   , char* buffer);
int   itoa_len(int i   , char* buffer);
void  dtoa    (double x, char* buffer, char* format, int numDecimal); // supported formats are "e" or "E" for scientific notation
int   dtoa_len(double x, char* buffer, char* format, int numDecimal); // or "f" or "F" for standard decimal notation.
void  stoa    (float x, char* buffer, char* format, int numDecimal); // in either case, numDecimal specifies number of chars after
int   stoa_len(float x, char* buffer, char* format, int numDecimal); // the decimal point
