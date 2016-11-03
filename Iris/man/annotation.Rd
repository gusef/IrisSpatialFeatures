\name{annotation}
\alias{annotation}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("annotation")}
\format{
  A data frame with 37 observations on the following 39 variables.
  \describe{
    \item{\code{SUBJID}}{a numeric vector}
    \item{\code{patid}}{a factor with levels \code{ML1403226} \code{ML1403234} \code{ML1403235} \code{ML1403236} \code{ML1403237} \code{ML1403238} \code{ML1403239} \code{ML1403240} \code{ML1403241} \code{ML1403242} \code{ML1403255} \code{ML1403261} \code{ML1403262} \code{ML1403263} \code{ML1403266} \code{ML1403267} \code{ML1403268} \code{ML1403269} \code{ML1403271} \code{ML1403272} \code{ML1403273} \code{ML1403274} \code{ML1403275} \code{ML1403276} \code{ML1403277} \code{ML1403280} \code{ML1403281} \code{ML1403283} \code{ML1408233} \code{ML1408236} \code{ML1408237} \code{ML1408239} \code{ML1408244} \code{ML1408245} \code{ML1408246} \code{ML1408248} \code{ML1408249} \code{ML1408250} \code{ML1408251} \code{ML1408253} \code{ML1408254} \code{ML1408255} \code{ML1408256} \code{ML1408259} \code{ML1408260} \code{ML1408261} \code{ML1408262} \code{ML1408264} \code{ML1408265} \code{ML1408267} \code{ML1408269} \code{ML1408270} \code{ML1408272} \code{ML1408273} \code{ML1408274} \code{ML1408276} \code{ML1408277} \code{ML1408278} \code{ML1408279} \code{ML1408280} \code{ML1408281} \code{ML1408282} \code{ML1408283} \code{ML1408284} \code{ML1408288} \code{ML1409107} \code{ML1409108} \code{ML1409110} \code{ML1409113} \code{ML1409122} \code{ML1409126} \code{ML1409127} \code{ML1409128} \code{ML1504117} \code{ML1504118} \code{ML1504119} \code{ML1504120} \code{ML1504121} \code{ML1504122} \code{ML1504123} \code{ML1504124} \code{ML1504126} \code{ML1504129} \code{ML1504131} \code{ML1504133} \code{ML1504134} \code{ML1504135} \code{ML1504137} \code{ML1504139} \code{ML1504140} \code{ML1504141} \code{ML1504144} \code{ML1504146} \code{ML1504147} \code{ML1504148} \code{ML1504150} \code{ML1504151} \code{ML1504152} \code{ML1504153} \code{ML1504155} \code{ML1504157} \code{ML1504158} \code{ML1504159} \code{ML1504160} \code{ML1504161} \code{ML1504162} \code{ML1504163} \code{ML1504164} \code{ML1504165} \code{ML1504166} \code{ML1504167} \code{ML1504168} \code{ML1504169}}
    \item{\code{b2m_IHC_p}}{a numeric vector}
    \item{\code{b2m_int}}{a numeric vector}
    \item{\code{b2m_hscore}}{a numeric vector}
    \item{\code{b2m_cat}}{a numeric vector}
    \item{\code{mhc1_ihc_p}}{a numeric vector}
    \item{\code{mhc1_int}}{a numeric vector}
    \item{\code{mhc1_hscore}}{a numeric vector}
    \item{\code{mhc1_cat}}{a numeric vector}
    \item{\code{mhc2_ihc_p}}{a numeric vector}
    \item{\code{mhc2_int}}{a numeric vector}
    \item{\code{mhc2_hscore}}{a numeric vector}
    \item{\code{mhc2_cat}}{a numeric vector}
    \item{\code{SITEID}}{a numeric vector}
    \item{\code{TRTP}}{a factor with levels \code{IPILIMUMAB-NIVOLUMAB} \code{NIVOLUMAB-IPILIMUMAB}}
    \item{\code{TRTPN}}{a numeric vector}
    \item{\code{ACTARM}}{a factor with levels \code{IPILIMUMAB-NIVOLUMAB} \code{NIVOLUMAB-IPILIMUMAB}}
    \item{\code{ACTARMCD}}{a factor with levels \code{A} \code{B}}
    \item{\code{RANDFL}}{a factor with levels \code{Y}}
    \item{\code{SAFFL}}{a factor with levels \code{Y}}
    \item{\code{RANDDT}}{a Date}
    \item{\code{TR01SDT}}{a Date}
    \item{\code{RESW13}}{a factor with levels \code{PD} \code{PR} \code{SD}}
    \item{\code{PROGW13}}{a factor with levels \code{N} \code{Y}}
    \item{\code{RESW13I}}{a factor with levels \code{PD} \code{PR} \code{SD}}
    \item{\code{RESW25}}{a factor with levels \code{CR} \code{CRU} \code{PD} \code{PR} \code{PRU} \code{SD}}
    \item{\code{PROGW25}}{a factor with levels \code{N} \code{Y}}
    \item{\code{BOR}}{a factor with levels \code{CR} \code{NE} \code{PD} \code{PR} \code{SD}}
    \item{\code{BORI}}{a factor with levels \code{CR} \code{PD} \code{PR} \code{SD}}
    \item{\code{FSTPDDT}}{a Date}
    \item{\code{DTHDT}}{a Date}
    \item{\code{DTHDTC}}{a Date}
    \item{\code{PFSFL41}}{a factor with levels \code{Y}}
    \item{\code{PFSFL71}}{a factor with levels \code{Y}}
    \item{\code{ADT}}{a Date}
    \item{\code{os_cns}}{a numeric vector}
    \item{\code{os_mos}}{a numeric vector}
    \item{\code{mhc2c}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(annotation)
## maybe str(annotation) ; plot(annotation) ...
}
\keyword{datasets}
