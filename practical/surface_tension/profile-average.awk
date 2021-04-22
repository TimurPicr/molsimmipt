NR < 3 {next}
NR == 4 {
  nitems = $2
  for (i = 1; i<=nitems; i++)
    n[i] = 0
    vx[i] = 0
  nlskip = (nitems + 1) * nskip + 3
}
NR <= nlskip {next}
{
  if (($1 != "#") && ($1 <= nitems))
    n[$1] += 1
    output[$1] += $ncol
}
END {
  for (i = 1; i<=nitems; i++)
  printf("%d %d %g\n",i, n[i], output[i]/n[i])
}
