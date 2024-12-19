source("general/general.R")





source("general/GenericFigures.R")

################################################################################
# FIGURE 1
################################################################################

# 1C____________________________________________________________________________
PlotEquilSiteWithInput <- function(mirna, column, experiment="equilibrium",
                                   n_constant=5, sitelist="resubmissionfinal",
                                   uniq=FALSE, combined=TRUE, singleonly=TRUE,
                                   height=4.5, width=4.5, buffer=FALSE,
                                   pdf.plot=FALSE) {
  sXc <- SubfunctionCall(SitesXCounts)
  x <- Norm(sXc[,1 + combined])
  y <- Norm(sXc[,column])
  site.colors <- kSiteColors[rownames(sXc)]
  SubfunctionCall(FigureSaveFile)
  xmin <- 5e-4
  xmax <- 1
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log="xy")
  xmax <- 1.00
  # Make x=y line:
  segments(xmin, ymin, xmax, ymax, lty=line_dash_length)
  # Make the lines connecting the points to the x = y line:
  segments(x, x, x, y, lty=line_dash_length, col=site.colors)
  # Make axes:
  AddLogAxis(1, label="Input library (%)", percent=TRUE)
  AddLogAxis(2, label="AGO-bound library (%)", percent=TRUE)
  # Add the points to the plot:
  Points(x, y, col=site.colors)
  R <- y/x
  names(R) <- rownames(sXc)
  A.stock <- kAgoStock[mirna, "equilibrium"] 
  A <- as.numeric(colnames(sXc)[column])/100*A.stock
  A.pM <- round(A*1000, 1)
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy')
  ago_text <- " pM AGO2-"
  text(xy[1], xy[2], labels=paste0(A.pM,ago_text, mirna),
       adj=c(0, 1))
  legend.coords <- GetPlotFractionalCoords(0.95, 0.025, log='xy')

  Legend(legend.coords, legend=ConvertTtoUandMmtoX(rownames(sXc)),
         col=kSiteColors[rownames(sXc)], xjust=1, yjust=0, y.intersp=0.9)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}





## FIGURES FOR RBNS EQUILIBRIUM PAPER ##########################################
MakeFigure1 <- function(uniq=FALSE) {
  message("Making Fig. 1")
  PlotEquilSiteWithInput("miR-1", 7, sitelist="canonical", combined=FALSE,
                         buffer=TRUE, pdf.plot="1.C")
  message("Done 1.C")
  # PlotSiteEnrichments("miR-1", sitelist="canonical", combined=FALSE,
  #                     buffer=TRUE, write_kds=TRUE, pdf.plot="1.D")
  # message("Done 1.D")
  # PlotSiteEnrichments("miR-1", combined=FALSE, remove_sites=FALSE, buffer=TRUE,
  #                     pdf.plot="1.E")
  # message("Done 1.E")
  # PlotSiteKds("miR-1", combined=FALSE, remove_sites=FALSE,
  #             buffer=TRUE, pdf.plot="1.F")
  # message("Done 1.F")
  # PlotSiteOccupancy("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="1.G")
  # message("Done Fig. 1")
}

MakeFigure2 <- function() {
  message("Making Fig. 2")
  PlotSiteKds("let-7a", adjusted_height=TRUE, pdf.plot="2.Ai")
  PlotSiteOccupancy("let-7a", adjusted_height=TRUE, pdf.plot="2.Aii")
  PlotSiteKds("miR-155", adjusted_height=TRUE, pdf.plot="2.Bi")
  PlotSiteOccupancy("miR-155", adjusted_height=TRUE, pdf.plot="2.Bii")
  PlotSiteKds("miR-124", collapse_AA=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ci")
  PlotSiteKds("miR-124", collapse_AA=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ci")
  PlotSiteOccupancy("miR-124", collapse_AA=TRUE, compcorrect=FALSE,
                    adjusted_height=TRUE, pdf.plot="2.Cii")

  message("Done Fig. 2")
}

MakeFigure3 <- function() {
  message("Making Fig. 3")
  PlotPositionalKds(pdf.plot="3.A")
  Plot8merVs7merCombined(pdf.plot="3.B")
  PlotSiteKdsVsSPS(pdf.plot="3.C")
  PlotSiteKdsVsRepression("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="3.D")
  PlotSiteKdsVsRepression("let-7a", pdf.plot="3.E")
  PlotSiteKdsVsRepression("miR-155", pdf.plot="3.F")
  PlotSiteKdsVsRepression("miR-124", pdf.plot="3.G")
  PlotSiteKdsVsRepression("lsy-6", pdf.plot="3.H")
  PlotSiteKdsVsRepression("miR-7-23nt", experiment="equilibrium2_nb",
                          pdf.plot="3.I")
  message("Done Fig. 3")
}

MakeFigure4 <- function() {
  message("Making Fig. 4")
  PlotSiteFlankEnrichments("miR-1", "8mer", combined=TRUE, combined_site=FALSE,
                           buffer=TRUE, pdf.plot="4.A")
  PlotSiteFlankKds("miR-1", combined=TRUE, combined_site=FALSE, buffer=TRUE,
                   pdf.plot="4.B")
  PlotFlankLinModel(pdf.plot="4.C_left")
  PlotFlankLinModelCoefficients(pdf.plot="4.C_right")
  PlotStructureVsFlankingKds("miR-1", "8mer", combined=TRUE, buffer=TRUE,
                             pdf.plot="4.D")
  message("Done Fig. 4")
}

MakeSupplementaryFigure1 <- function() {
  message("Making fig. S1")
  PlotAgoPrepPurity(pdf.plot="S1.A", no_marker=FALSE, no_adapter=FALSE)
  PlotMiR1_KmersCorrelation(pdf.plot="S1.B", kmer_len=9)
  PlotEnrichmentsAgainstKds("miR-1", sitelist="canonical", buffer=TRUE,
                            combined=FALSE, pdf.plot="S1.C")
  PlotPositionalEnrichment("miR-1", buffer=TRUE, pdf.plot="S1.Di") # checked
  PlotPositionalEnrichment("let-7a", pdf.plot="S1.Dii") # checked
  PlotPositionalEnrichment("miR-155", pdf.plot="S1.Diii") # checked
  PlotPositionalEnrichment("miR-124", pdf.plot="S1.Div") # checked
  PlotPositionalEnrichment("lsy-6", pdf.plot="S1.Dv") # checked
  PlotPositionalEnrichment("miR-7-23nt", exp="equilibrium2_nb", pdf.plot="S1.Dvi")
  PlotPositionalEnrichment("miR-155", sites=k3PSites, pdf.plot="S1.Ei") # checked
  PlotPositionalEnrichment("miR-124", sites=k3PSites, pdf.plot="S1.Eii") # checked
  PlotPositionalEnrichment("lsy-6", sites=k3PSites, pdf.plot="S1.Eiii") # checked  break
  PlotWorstSiteKdCrossValScatter("miR-1", sitelist="canonical", combined=FALSE,
                                 buffer=TRUE, pdf.plot="S1.F") # checked
  PlotSiteKdCrossValMatrix("miR-1", sitelist="canonical", combined=FALSE,
                           buffer=TRUE, pdf.plot="S1.G") # checked
  PlotSalomonComparison(pdf.plot="S1.H") # checked
  message("Done fig. S1")
}

MakeSupplementaryFigure2 <- function() {
  message("Making fig. S2.")
  PlotCompetitorOligoSiteSchematic(pdf.plot="S2.A")
  PlotSiteKds("lsy-6", adjusted_height=TRUE, pdf.plot="S2.Bi")
  PlotSiteOccupancy("lsy-6", sitelist="resubmissionfinal",  adjusted_height=TRUE,
                    pdf.plot="S2.Bii")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", adjusted_height=TRUE,
              pdf.plot="S2.Ci")
  PlotSiteOccupancy("miR-7-23nt", exp="equilibrium2_nb", adjusted_height=TRUE,
                    pdf.plot="S2.Cii")
  message("Done fig. S2.")
}

MakeSupplementaryFigure3 <- function() {
  message("Making fig. S3.")
  PlotBaekKds("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="S3.A")
  PlotBaekKds("let-7a", pdf.plot="S3.B")
  PlotBaekKds("miR-155", pdf.plot="S3.C")
  PlotBaekKds("miR-124", pdf.plot="S3.D")
  PlotBaekKds("lsy-6", pdf.plot="S3.E")
  PlotBaekKds("miR-7-23nt", experiment="equilibrium2_nb", pdf.plot="S3.F")
  message("Done fig. S3.")
}

MakeSupplementaryFigure4 <- function() {
  message("Making fig. S4.")
  PlotBulgeKds("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="S4.B")
  PlotDelKds("miR-1", , combined=FALSE, buffer=TRUE, pdf.plot="S4.C")
  PlotBulgeKds("let-7a", pdf.plot="S4.D")
  PlotDelKds("let-7a", pdf.plot="S4.E")
  PlotBulgeKds("miR-155", pdf.plot="S4.F")
  PlotDelKds("miR-155", pdf.plot="S4.G")
  PlotBulgeKds("miR-124", pdf.plot="S4.H")
  PlotDelKds("miR-124", pdf.plot="S4.I")
  PlotBulgeKds("lsy-6", pdf.plot="S4.J")
  PlotDelKds("lsy-6", pdf.plot="S4.K")
  PlotBulgeKds("miR-7-23nt", experiment="equilibrium2_nb", pdf.plot="S4.L")
  PlotDelKds("miR-7-23nt", experiment="equilibrium2_nb",  pdf.plot="S4.M")
  message("Done fig. S4.")
}


MakeSupplementaryFigure5 <- function() {
  message("Making fig. S5.")
  PlotReporterAssaySchematicLetters(pdf.plot="S5.A")
  PlotReporterAssayKdVsL2fc("miR-1", pdf.plot="S5.Bi")
  PlotReporterAssayKdVsL2fc("let-7a", pdf.plot="S5.Bii")
  PlotReporterAssayKdVsL2fc("miR-155", pdf.plot="S5.Biii")
  PlotReporterAssayKdVsL2fc("miR-124", pdf.plot="S5.Biv")
  PlotReporterAssayKdVsL2fc("lsy-6", pdf.plot="S5.Bv")
  PlotReporterAssayKdVsL2fc("miR-7", pdf.plot="S5.Bvi")
  message("Done fig. S5.")
 }

MakeSupplementaryFigure6 <- function() {
  message("Making fig. S6.")
  PlotSiteFlankKds("let-7a", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.A")
  PlotSiteFlankKds("miR-155", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.B")
  PlotSiteFlankKds("miR-124", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.C")
  PlotSiteFlankKds("lsy-6", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.D")
  PlotSiteFlankKds("miR-7-23nt", experiment="equilibrium2_nb", width=7.3,
                   adjusted_height=TRUE, pdf.plot="S6.E")
  PlotFlankKdsVsRepression(sitelist="topcanonical", pdf.plot="S6.F")
  PlotAllSiteInputDistribution("miR-1", "8mer", combined=TRUE, buffer=TRUE,
                               pdf.plot="S6.G")
  PlotFlankPlCorrelationWindowsNew("miR-1", "8mer",
                                   buffer=TRUE, pdf.plot="S6.H")
  PlotFlanksSamplePl("miR-1", "8mer", "0.4", buffer=TRUE, matchdist=TRUE,
                     pdf.plot="S6.I")
  message("Done fig. S6.")
 }

MakeRefereeResponseFigure <- function() {
  message("Making referee response fig.")
  PlotFlowthrough("let-7a", pdf.plot="R.A")
  PlotFlowthrough("let-7a", pdf.plot="R.B", log_y=FALSE)
  PlotFlowthrough("let-7a", pdf.plot="R.C",
                  enrich_against_enrich=TRUE)
  message("Done referee response fig.")
 }



 
MakeFigure1()
# MakeFigure2()
# MakeFigure3()
# MakeFigure4()
# MakeSupplementaryFigure1()
# MakeSupplementaryFigure2()
# MakeSupplementaryFigure3()
# MakeSupplementaryFigure4()
# MakeSupplementaryFigure5()
# MakeSupplementaryFigure6()



