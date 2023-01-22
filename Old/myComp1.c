int cmpfunc (const void * a, const void * b) {
   if ( *(float*)a < *(float*)b ) return -1;
   if ( *(float*)a > *(float*)b ) return 1;
   return(0);
}
