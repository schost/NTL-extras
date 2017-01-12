{
  if ($0 == ""){
    new = 1;
  }

  if ($3 > $4 && new == 1){
    print $0;
    new = 0;
  }
    
}
