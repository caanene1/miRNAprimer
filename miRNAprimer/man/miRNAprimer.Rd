\name{miRNAprimer}
\alias{miRNAprimer}
\title{qPCR primers for miRNA quantification
}
\description{ Derives forward primers and the unique part of the stemloop primer.
}
\usage{
miRNAprimer(x, s, r)
}
\arguments{
  \item{x}{Mature miRNA sequence.
}
  \item{s}{Base stemloop sequence.
}
  \item{r}{Universal reverse sequence.
}
}
\details{Stemloop RT qPCR is the gold standard for the quantification of individual miRNAs in RNA samples. The manual design of the three sets of primers can be challenging in large experiments. Many tools have been presented, but they often offer highly restrictive web tools and not available on R.
}

\value{Table of stemloop primer, reverse primer and several miRNA foward primers. The melting temperatures are included for primer selection.
}
\references{PMID: 21732315,
            PMID: 23383094
}
\author{Chinedu A. Anene
}
\note{The user need to obtain and provide the mature miRNA sequence in 5' to 3' orientation. Sequences in the examples from PMID: 23383094 can used. The next update will implement the download of mature miRNA seqeuance from micro.org. Other published stemloop will be listed in the next update. Example stemloop and reverse sequence are from PMID: 23383094.
}
\seealso{
}
\examples{
r <- "GTGCAGGGTCCGAGGT"
s <- "GTTGGCTCTGGTGCAGGGTCCGAGGTATTCGCACCAGAGCCAAC"
x <- "CGCAAGGATGACACGCAAATTC"
miRNAprimer(x,s,r)
}
\keyword{ ~primer

}
\keyword{ ~miRNA

}
