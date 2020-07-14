//https://groups.google.com/forum/#!topic/stan-users/94BQm3MKIqw
functions {

  matrix subset_stan(matrix x, int[] is, int target) { 
   int count; 
   if (size(is) != rows(x)) 
     reject("illegal input");  // stricter than R's broadcasting 

   // count number of rows in result 
   count = 0; 
   for (n in 1:size(is)) 
     if (is[n] == target) 
       count = count + 1; 

   // assign rows in result 
   { 
     matrix[count, cols(x)] result; 
     int pos; 
     pos = 1; 
     for (n in 1:size(is)) { 
       if (is[n] == target) { 
         result[pos] = x[n];  // assigns entire row 
         pos = pos + 1; 
       } 
     } 
     return result; 
   } 
  }
}
model {
}
