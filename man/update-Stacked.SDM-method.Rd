% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.stack.R
\name{update,Stacked.SDM-method}
\alias{update,Stacked.SDM-method}
\title{Update a previous SSDM}
\usage{
\S4method{update}{Stacked.SDM}(
  object,
  Occurrences,
  Env,
  Xcol = "Longitude",
  Ycol = "Latitude",
  Pcol = NULL,
  Spname = NULL,
  name = stack@name,
  save = FALSE,
  path = getwd(),
  thresh = 1001,
  tmp = FALSE,
  verbose = TRUE,
  GUI = FALSE,
  ...
)
}
\arguments{
\item{object}{Stacked.SDM. The previously built SSDM.}

\item{Occurrences}{data frame. New or updated occurrence table (can be
processed first by \code{\link{load_occ}}).}

\item{Env}{raster object. Environment raster object (can be processed first by
\code{\link{load_var}}).}

\item{Xcol}{character. Name of the column  in the occurrence table  containing
Latitude or X coordinates.}

\item{Ycol}{character. Name of the column in the occurrence table  containing
Longitude or Y coordinates.}

\item{Pcol}{character. Name of the column in the occurrence table specifying
whether a line is a presence or an absence. A value of 1 is presence and
value of 0 is absence. If NULL presence-only dataset is assumed.}

\item{Spname}{character. Name of the new or updated species.}

\item{name}{character. Optional name given to the final SSDM produced, by
default it's the name of the previous SSDM.}

\item{save}{logical. If set to true, the model is automatically saved.}

\item{path}{character. Name of the path to the directory to contain the saved
SSDM.}

\item{thresh}{numeric. A single integer value representing the number of equal
interval threshold values between 0 and 1.}

\item{tmp}{logical. If set to true, the habitat suitability map of each
algorithm is saved in a temporary file to release memory. But beware: if you
close R, temporary files will be deleted To avoid any loss you can save
your model with \code{\link{save.model}}.}

\item{verbose}{logical. If set to true, allows the function to print text in
the console.}

\item{GUI}{logical. Don't take that argument into account (parameter for the
user interface).}

\item{...}{additional parameters for the algorithm modelling function (see
details below).}
}
\value{
an S4 \linkS4class{Stacked.SDM} class object viewable with the
 \code{\link{plot.model}} function.
}
\description{
Update a previous SSDM with new occurrence data. The function takes as inputs
updated or new occurrence data from one species, previous environmental
variables, and an S4 \linkS4class{Stacked.SDM} class object containing a
previously built SSDM.
}
\examples{
\dontrun{
update(stack, Occurrences, Env, Spname = 'NewSpecie')
}

}
\seealso{
\code{\link{stack_modelling}} to build SSDMs.
}
