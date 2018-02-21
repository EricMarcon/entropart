ArgumentOriginalName <- 
function(x) 
{
  # Code from https://stackoverflow.com/questions/26557731/, modified.
  
  # First get the name of the argument received as "x"
  Received <- quote(substitute(x))
  VarName <- eval(Received)

  # Search the frames up to the top one to find the real name of the variable or expression
  for(i in rev(utils::head(sys.frames(), -1L))) {
    Received[[2]] <- VarName
    # Candidate name
    VarNameCandidate <- eval(Received, i)
    VarNameChar <- as.character(VarNameCandidate)
    if(length(VarNameChar)) {
      # VarNameCandidate must not be NULL
      if(VarNameChar[length(VarNameChar)] != "NULL") {
        # If Received is a dataframe, then the top frame may return Name$NULL. Else, validate the candidate.
        VarName <- VarNameCandidate
      }
    }
  }
  argumentOriginalName <- deparse(VarName)
  # Not found or too long
  if (length(argumentOriginalName) > 1) {
    argumentOriginalName <- "-"
  } else {
    if (nchar(argumentOriginalName) > 30)
      paste(substr(argumentOriginalName, 1, 25), "(...)", sep="")
  }
  return(argumentOriginalName) 
}
