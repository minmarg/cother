%module ddms
%{
    #include "ddms.h"  
%}

%include "carrays.i"
%array_functions(int, intArray);

%include "ddms.h"
