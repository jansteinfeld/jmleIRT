# .onAttach <- function(libname, pkgname) {
#   supports_color <- function() {
#     # Conservative check; base R only
#     nzchar(Sys.getenv("TERM")) && !identical(getOption("crayon.enabled", FALSE), FALSE)
#   }

#   clr <- list(
#     accent = "\033[38;5;39m",   # blue
#     faint  = "\033[2m",         # faint
#     reset  = "\033[0m"
#   )
#   C <- if (supports_color()) clr else list(accent = "", faint = "", reset = "")

#   ver  <- tryCatch(utils::packageVersion(pkgname), error = function(e) "unknown")
#   rver <- paste0(R.version$major, ".", R.version$minor)

#   title <- sprintf("%s%s%s %s", C$accent, pkgname, C$reset, ver)
#   line  <- paste0(C$faint, strrep("-", 62), C$reset)

#   body <- paste0(
#     "R ", rver, " | ?",
#     pkgname, " | citation('", pkgname, "')"
#   )

#   msg <- paste0(
#     "\n", line, "\n",
#     " ", title, "\n",
#     " ", body, "\n",
#     " Rasch JML (C++). optional Warm's WLE  NA-tolerant.  centered\n",
#     line, "\n"
#   )
#   packageStartupMessage(msg)
# }


# in R/zzz.R
.onAttach <- function(libname, pkgname) {
  pkg  <- pkgname
  ver  <- tryCatch(utils::packageVersion(pkg), error = function(e) "unknown")
  rver <- paste0(R.version$major, ".", R.version$minor)

  msg <- c(
    "",
    "  |--------------------------------------------------------------|",
    sprintf("  |  %s %s   |", format(sprintf("%s %s", pkg, ver), width = 56, justify = "left"), strrep(" ", 0)),
    sprintf("  |  R %s | type ?%s for help | citation(' %s ')   |",
            format(rver, width = 6, justify = "right"),
            pkg,
            pkg),
    "  |                                                              |",
    "  |  Rasch (1PL) JML with optional Warm's WLE and NA handling.   |",
    "  |  Fast C++ core, centering and epsilon adjustment.            |",
    "  |--------------------------------------------------------------|",
    ""
  )
  packageStartupMessage(paste(msg, collapse = "\n"))
}
