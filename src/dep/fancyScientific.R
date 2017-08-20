
# ###
# Fancy scientific notation for axes, taken and modified from Brian Diggs 
# https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4
# ###

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE, digits = 2)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}