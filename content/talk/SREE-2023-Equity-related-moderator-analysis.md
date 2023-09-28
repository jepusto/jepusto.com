---
publishDate: 2023-09-27T16:00:00  # Schedule page publish date.
title: "Equity-related moderator analysis in syntheses of dependent effect sizes: Conceptual and statistical considerations"
date: 2023-09-27T16:15:00
event: "Society for Research on Educational Effectiveness Conference"
event_short: "SREE 2023"
location: "Arlington, VA"
featured: false
url_slides: "/files/SREE-2023-Equity-Related-Moderator-Analysis.pdf"
math: false
highlight: false
abstract: ""
event_url: "https://www.sree.org/2023-conference"
authors: [admin, Jingru Zhang, Elizabeth Tipton]
projects: 
summary: 
tags: []
---

# Background/Context

In meta-analyses examining educational interventions, researchers seek to understand the distribution of intervention impacts, in order to draw generalizations about what works, for whom, and under what conditions. One common way to examine equity implications in such reviews is through moderator analysis, which involves modeling how intervention effect sizes vary depending on the characteristics of primary study participants. For example, one might estimate associations between effect size and the percentage of the primary study participants who were from a rural school, from a low-income family, identified as a specific racial or ethnic group, or designated as an English Language Learner. Such moderator analyses can provide insights about the populations and contexts where an intervention is more or less effective—that is, they can address questions of who benefits and how the effects of an educational intervention are distributed. 

Meta-analyses of educational interventions often include primary studies that report multiple relevant effect size estimates, such for more than one measure of an outcome construct, at multiple time-points, for multiple versions of an intervention, or for different sub-groups of participants. This leads to a data structure where the effects from a given study are correlated, necessitating the use of statistical methods that are appropriate for dependent observations. Methodological research in this area has provided estimation and inference methods that which can handle dependent effect sizes, including multi-level meta-analyses (Van den Noortgate et al., 2013, 2015), robust variance estimation (Hedges et al., 2010), and combinations thereof (Fernández-Castilla et al., 2020; Pustejovsky & Tipton, 2022). However, there has been much less attention to the specific forms of moderator analysis that are of interest in practice. 

# Purpose/Objective/Research Question

We aim to identify conceptual and statistical considerations for moderator analysis of equity-related variables in meta-analyses involving dependent effect sizes. Specifically, we distinguish between direct evidence and contextual evidence about equity of impacts and show that the choice of meta-analytic model can be consequential for analyses involving direct evidence. We then examine how meta-analysts currently conduct equity-related moderator analyses, by reviewing completed research synthesis projects funded by the Institute of Education Sciences (IES) over the period of 2002 to 2018. We find that most projects do not distinguish between direct and contextual evidence and use analytic approaches that are inefficient for synthesizing direct evidence.
Conceptual Considerations

Moderator analyses of equity-related variables can be carried out by regressing effect size estimates on predictors encoding participant characteristics. Consider a synthesis in which some primary studies contribute multiple effect sizes. In this data structure, a predictor might represent a study-level characteristic or one that varies across the effects within a given study. The level of variation is especially salient for analysis of equity-related variables because study-level characteristics and effect-level characteristics represent qualitatively different types of evidence. For study-level predictors, associations with effect size pertain to the study’s context and are not necessarily indicative of individual-level variation in impacts. Thus, interpretation is challenging because studies vary in many ways, with many possible sources of confounding. For effect-level predictors, within-study variation represents direct evidence about individual-level moderation (e.g., a comparison of impacts between low-income and higher-income participants in the same study), unconfounded by study-level characteristics.

We describe different strategies for separately investing direct and contextual evidence about moderation, including a) decomposing the predictor into study-level average and within-study centered components or b) inclusion of the raw predictor and the study-level average in a meta-regression. Although strategy (a) has been recommended previously in the context of meta-analysis of dependent effects (Tanner-Smith & Tipton, 2014), our presentation makes explicit the connection to equity-related moderator analysis and specifies the data requirements for applying it. 

# Statistical Considerations

Meta-regression with dependent effect sizes involves choosing a working model for the dependence structure, which determines the set of weights used for estimating the meta-regression. Several different working models have been proposed, including correlated effects and hierarchical effects models (Hedges et al., 2010), a correlated-and-hierarchical effects model (Pustejovsky & Tipton, 2022) and the multi-level meta-analysis model (Van den Noortgate et al., 2013, 2015). Ad hoc strategies, such as aggregating effects to the study level or ignoring dependence, can also be understood as working models. 

Previous research and guidance about the choice of working model has argued that the choice of working model is fairly inconsequential so long as the working model is roughly similar the true dependence structure (Hedges et al., 2010; Tanner-Smith et al., 2016; Tanner-Smith & Tipton, 2014). In the appendix, we examine the exact weights assigned by a variety of different working models to studies with direct evidence and contextual evidence. Contrary to past guidance, we find that different working models can lead to quite different weighting—particularly for direct evidence (i.e., study mean-centered predictors).

# Current Practice

To understand current practices for analysis of equity-related moderator variables, we reviewed completed meta-analysis projects funded by IES over the period of 2002 to 2018. We identified grants that (a) had journal articles reporting a meta-analysis, (b) were not methodological, and (c) were not training programs. A search of the IES website for project descriptions that included the word “meta-analysis” returned 80 results, of which 25 met inclusion criteria. Table 1 summarizes the approaches to moderator analyses used in these projects. Most projects reported some form of meta-regression analysis, but very few described a centering strategy and only one project used study-mean centering. Notably, the correlated effects and hierarchical effects working models were commonly used, yet these models involve inefficient weighting of direct evidence.  

# Conclusions

In light of this review of current practice, the conceptual and statistical considerations that we describe suggest that there is substantial room for improvement in how meta-analysts conduct moderation analysis, particularly for equity-related variables where individual-level variation is of primary interest. Even under this simple—simplistic, even—conception of equity, bringing systematic review and meta-analysis methods to bear to address inequities in the education system will require not only improving analytic practices, but also changing how primary investigations frame questions, collect data, and report findings.