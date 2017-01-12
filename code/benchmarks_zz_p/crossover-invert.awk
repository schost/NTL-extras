BEGIN{
  go = 1;
  inc = 0;
}

{
  if (inc == 1 && NF > 0){
    if ($3 <= t0){
      n = $1;
      a = $2;
    }
  }
  if (go == 1){
    t0 = $3;
    go = 0;
    inc = 1;
  }
  if (NF == 0){
    go = 1;
    printf("%i %i\n", n, a);
    inc = 0;
  }
}
