GenSimpson <-
function(Ps, r, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (sum(Ps*(1-Ps)^r))
}
