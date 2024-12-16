source("general/general.R")
source("general/general_kinetics.R")
library(wCorr)
graphics.off()


pt_cex_final <- 1.35
# pch21_lwd.cex <- 0.7 * par()$cex
legend_pch <- 19
line_dash_length <- 2
occupancy_schematic_bg_cex <- 0.45
occupancy_schematic_bg_col <- "#989898"
occupancy_mirna_seq_col <- "#958872"

AGO_mir_label <- "[AGO2-miRNA] (pM)"
AGO_mir1_label <-"[AGO2-miR-1]"
AGO_mir1_label_no_conc <- "AGO2-miR-1 library"





FigureSaveFile <- function(pdf.plot, height=5, width=5, xpos=20, ypos=20) {
  if (class(pdf.plot) == "character") {
    print(pdf.plot)
    figure.split <- unlist(strsplit(pdf.plot, split="\\.", perl=TRUE))
    figure <- figure.split[1]
    if (length(figure.split) > 2) {
      panel <- paste0(figure.split[2:length(figure.split)], collapse=".")
    } else {
      panel <- figure.split[2]    
    }
    figurename <- paste0(figure, panel)
    path <- paste0("KineticsPaper/Figure_", figure, "/", figurename, "_raw.pdf")
    pdf(file=path, height=2.2/5*height, width=2.2/5*width, useDingbats=FALSE)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }  
}


test <- GetKineticsData("miR-1")


break

################################################################################
# FIGURE 1
################################################################################

# 1C____________________________________________________________________________

## FIGURES FOR RBNS KINETICS PAPER ##########################################

MakeFigure1 <- function(uniq=FALSE) {
  message("Making Fig. 1")
  PlotEquilSiteWithInput("miR-1", 7, sitelist="canonical", combined=FALSE,
                         buffer=TRUE, pdf.plot="1.C")
  PlotSiteEnrichments("miR-1", sitelist="canonical", combined=FALSE,
                      buffer=TRUE, pdf.plot="1.D")
  PlotSiteEnrichments("miR-1", combined=FALSE,  buffer=TRUE, pdf.plot="1.E")
  PlotSiteKds("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="1.F")
  PlotSiteOccupancy("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="1.G")
  message("Done Fig. 1")
}

MakeFigure2 <- function() {
  message("Making Fig. 2")
  PlotSiteKds("let-7a", sitelist="paperfinal", trim=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ai")
  PlotSiteOccupancy("let-7a", sitelist="paperfinal", trim=TRUE,
                    adjusted_height=TRUE, pdf.plot="2.Aii")
  PlotSiteKds("miR-155", sitelist="paperfinal", trim=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Bi")
  PlotSiteOccupancy("miR-155", sitelist="paperfinal", trim=TRUE,
                    adjusted_height=TRUE, pdf.plot="2.Bii")
  PlotSiteKds("miR-124", sitelist="paperfinal", trim=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ci")
  PlotSiteOccupancy("miR-124", sitelist="paperfinal", trim=TRUE,
                    adjusted_height=TRUE, pdf.plot="2.Cii")
  PlotSiteKds("lsy-6", sitelist="paperfinal", trim=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Di")
  PlotSiteOccupancy("lsy-6", sitelist="paperfinal", trim=TRUE,
                    adjusted_height=TRUE, pdf.plot="2.Dii")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", combined=FALSE,
              sitelist="paperfinal", trim=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ei")
  PlotSiteOccupancy("miR-7-23nt", exp="equilibrium2_nb", combined=FALSE,
                    sitelist="paperfinal",
                    trim=TRUE, adjusted_height=TRUE,
                    pdf.plot="2.Eii")
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
                          combined=FALSE,
                          pdf.plot="3.I")
  message("Done Fig. 3")
}

MakeFigure4 <- function() {
  message("Making Fig. 4")
  PlotSiteFlankEnrichments("miR-1", "8mer", combined=TRUE, buffer=TRUE,
                           pdf.plot="4.A")
  PlotSiteFlankKds("miR-1", combined=TRUE, buffer=TRUE, pdf.plot="4.B")
  PlotFlankLinModel(pdf.plot="4.C_left")
  PlotFlankLinModelCoefficients(pdf.plot="4.C_right")
  PlotStructureVsFlankingKds("miR-1", "8mer", combined=TRUE, buffer=TRUE,
                             pdf.plot="4.D")
  PlotAllSamplePlFlanks("miR-1", "8mer", "0.4", buffer=TRUE, combined=TRUE,
                        pdf.plot="4.E")
  message("Done Fig. 4")
}

MakeSupplementaryFigure1 <- function() {
  message("Making fig. S1")
  PlotAgoPrepPurity(pdf.plot="S1.A", no_marker=FALSE, no_adapter=FALSE)
  PlotMiR1_KmersCorrelation(pdf.plot="S1.B", kmer_len=9)
  PlotEnrichmentsAgainstKds("miR-1", sitelist="canonical", buffer=TRUE,
                            combined=FALSE, pdf.plot="S1.C")
  PlotPositionalEnrichment("miR-1", buffer=TRUE, pdf.plot="S1.Di")
  PlotPositionalEnrichment("let-7a", pdf.plot="S1.Dii")
  PlotPositionalEnrichment("miR-155", pdf.plot="S1.Diii")
  PlotPositionalEnrichment("miR-124", pdf.plot="S1.Div")
  PlotPositionalEnrichment("lsy-6", pdf.plot="S1.Dv")
  PlotPositionalEnrichment("miR-7-23nt", exp="equilibrium2_nb", pdf.plot="S1.Dvi")
  PlotPositionalEnrichment("miR-155", sites=k3PSites, pdf.plot="S1.Ei")
  PlotPositionalEnrichment("miR-124", sites=k3PSites, pdf.plot="S1.Eii")
  PlotPositionalEnrichment("lsy-6", sites=k3PSites, pdf.plot="S1.Eiii")
  PlotWorstSiteKdCrossValScatter("miR-1", sitelist="canonical", combined=FALSE,
                                 buffer=TRUE, pdf.plot="S1.F")
  PlotSiteKdCrossValMatrix("miR-1", sitelist="canonical", combined=FALSE,
                           buffer=TRUE, pdf.plot="S1.G")
  PlotSalomonComparison(pdf.plot="S1.H")
  message("Done fig. S1")
}

MakeSupplementaryFigure2 <- function() {
  message("Making fig. S2.")
  PlotSiteKds("let-7a", adjusted_height=TRUE, pdf.plot="S2.A")
  PlotSiteKds("miR-155", adjusted_height=TRUE, pdf.plot="S2.B")
  PlotSiteKds("miR-124", adjusted_height=TRUE, pdf.plot="S2.C")
  PlotSiteKds("lsy-6", adjusted_height=TRUE, pdf.plot="S2.D")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", combined=FALSE,
              adjusted_height=TRUE, pdf.plot="S2.E")
  message("Done fig. S2.")
}

MakeSupplementaryFigure3 <- function() {
  message("Making fig. S3.")
  PlotBaekKds("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="S3.A")
  PlotBaekKds("let-7a", pdf.plot="S3.B")
  PlotBaekKds("miR-155", pdf.plot="S3.C")
  PlotBaekKds("miR-124", pdf.plot="S3.D")
  PlotBaekKds("lsy-6", pdf.plot="S3.E")
  PlotBaekKds("miR-7-23nt", combined=FALSE, experiment="equilibrium2_nb",
              pdf.plot="S3.F")
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
  PlotBulgeKds("miR-7-23nt", experiment="equilibrium2_nb", combined=FALSE,
               pdf.plot="S4.L")
  PlotDelKds("miR-7-23nt", experiment="equilibrium2_nb", combined=FALSE,
             pdf.plot="S4.M")
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
  message("Making fig. S5.")
  # PlotSiteFlankKds("let-7a", adjusted_height=TRUE, trim=TRUE, width=7.3,
  #                  pdf.plot="S6.A")
  # PlotSiteFlankKds("miR-155", adjusted_height=TRUE, trim=TRUE, width=7.3,
  #                  pdf.plot="S6.B")
  # PlotSiteFlankKds("miR-124", adjusted_height=TRUE, trim=TRUE, width=7.3,
  #                  pdf.plot="S6.C")
  # PlotSiteFlankKds("lsy-6", adjusted_height=TRUE, trim=TRUE, width=7.3,
  #                  pdf.plot="S6.D")
  # PlotSiteFlankKds("miR-7-23nt", experiment="equilibrium2_nb", width=7.3,
  #                  adjusted_height=TRUE, combined=FALSE, trim=TRUE,
  #                  pdf.plot="S6.E")
  # PlotFlankKdsVsRepression(sitelist="topcanonical", pdf.plot="S6.F")
  # PlotAllSiteInputDistribution("miR-1", "8mer", combined=TRUE, buffer=TRUE,
  #                              pdf.plot="S6.G")
  PlotFlankPlCorrelationWindowsNew("miR-1", "8mer",
                                   buffer=TRUE, pdf.plot="S6.H")
  # PlotFlanksSamplePl("miR-1", "8mer", "0.4", buffer=TRUE, matchdist=TRUE,
  #                    pdf.plot="S6.I")
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


