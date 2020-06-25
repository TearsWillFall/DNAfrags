#' Installs all the required tools to extract DNA fragment information.
#'
#' This function downloads and compiles the source files of all the
#' tools needed to manipulate DNA fragment data. However, it DOESN'T provide
#' all the libraries and dependencies necessary for their succesful installation,
#' therefore, they have to be installed beforehand as described in the README.md
#' or at https://github.com/TearsWillFall/DNAfrags.
#' @export


install_required_tools=function(){
  BTools::install_tools(whitelist=c("samtools"))
}
