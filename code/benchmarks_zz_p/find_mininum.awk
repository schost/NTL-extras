{
  if ($2 == 2000) {
    m = $4;
    idx = $3;
    for (i = 6; i < NF; i+=2){
      j = i-1;
      if ($i < m){
	m = $i;
	idx = $j;
      }
    }
    printf("%i %i %i\n", $1, $2, idx);
  }
}
