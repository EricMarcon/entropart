GenSimpsonD <-
function(Ps, r, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (1/(1-(sum(Ps*(1-Ps)^r))^(1/r)))
}
