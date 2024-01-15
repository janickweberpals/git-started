# README

## Background

This repository complements the `Let's git it started` introduction to
git for reproducible analytic workflows in real-world evidence (RWE)
studies.

<u>Note</u>: This is an example repository based on the HARPER
template<sup>1</sup> (available under
[gitlab-scm.partners.org/drugepi/harper](https://gitlab-scm.partners.org/drugepi/harper))
that uses a toy example to demonstrate the composition of a potential
RWE study repository and how to version-control and track changes.

## Mock study

### Dataset

For this mock study we will use the `smdi_data_complete` dataset that
comes with the smdi R package.<sup>2</sup>

### Background

``` r
library(smdi)
suppressPackageStartupMessages(library(dplyr))

# smdi_data_complete does not come with
# a patient ID, so we artificially create
# one for this example
data <- smdi_data_complete |> 
  mutate(patientID = as.character(paste0("patientID", dplyr::row_number())))

data |> 
  glimpse()
```

    Rows: 2,500
    Columns: 15
    $ exposure      <int> 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0,…
    $ age_num       <dbl> 35.24, 51.18, 88.17, 50.79, 40.52, 64.57, 73.58, 42.38, …
    $ female_cat    <int> 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1,…
    $ ecog_cat      <int> 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1,…
    $ smoking_cat   <int> 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1,…
    $ physical_cat  <int> 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,…
    $ egfr_cat      <int> 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1,…
    $ alk_cat       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    $ pdl1_num      <dbl> 45.03, 43.54, 41.74, 45.51, 31.28, 38.01, 47.28, 37.28, …
    $ histology_cat <int> 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,…
    $ ses_cat       <chr> "2_middle", "3_high", "2_middle", "2_middle", "2_middle"…
    $ copd_cat      <int> 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1,…
    $ eventtime     <dbl> 5.000000000, 4.754220474, 0.253391563, 5.000000000, 5.00…
    $ status        <int> 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1,…
    $ patientID     <chr> "patientID1", "patientID2", "patientID3", "patientID4", …

Let’s assume we want to conduct a RWE study comparing two antineoplastic
treatments in the context of real-world overall survival (rwOS)
outcomes.

### Eligibility criteria

For this example we want to apply a couple of eligibility criteria,
i.e.,

1.  Age \> 21 years at index date (`age_num > 21`)

2.  No history of smoking (`smoking_cat == 0`)

3.  PD-L1 biomarker expression \> 50% (`pdl1_num > 50`)

``` r
library(visR)

attrition <- get_attrition(
  # dataset
  data = data,
  
  # labels                               
  criteria_descriptions = c("Age > 21 years at index date",
                             "No history of smoking"),
  # coded condition
  criteria_conditions   = c("age_num > 21",
                             "smoking_cat == 0"),
  # patient ID column
  subject_column_name   = "patientID"
  
  )

visr(
  x = attrition,
  description_column_name = "Criteria", 
  value_column_name = "Remaining N"
  )
```

<img src="README_files/figure-commonmark/fig-attrition-1.png"
id="fig-attrition"
alt="Figure 1: Cohort attrition and final study size." />

### Propensity score matching

We perform a 1:1 propensity score matching<sup>3</sup>, adjusting for
pre-exposure covariates to control for confounding bias.

``` r
library(MatchIt)
suppressPackageStartupMessages(library(cobalt))
# specify ps model
exposure_form <- as.formula("exposure ~ age_num + female_cat + ses_cat")

# estimate 1:1 propensity score matching
# with 0.2 caliper on propensity score
ps_matching <- matchit(
  formula = exposure_form,
  data = data,
  ratio = 1,
  method = "nearest",
  distance = "glm",
  link = "logit",
  caliper = 0.2,
  replace = F
  )

# visualize covariate balance
love.plot(
  x = ps_matching,
  abs = TRUE, # if absolute values should be plotted
  thresholds = 0.1, # vertical line with "balance" threshold
  drop.distance = TRUE, # should the distance measure be removed
  var.order = "unadjusted", # variable order on y-axis
  colors = c("blue", "orange"), # first color =  unadjusted sample, second = adjusted sample
  stars = "std", # to indicate mean differences that have been standardized
  shapes = 17, # the shape of the SMD geometries
  size = 4, # the size of the SMD gemotries
  grid = TRUE,
  position = "top" # legend position
  )
```

![Covariate balance before and after propensity score
matching.](README_files/figure-commonmark/ps-matching-1.png)

### Estimating marginal hazard ratio

``` r
# use propensity score matched cohort
data_matched <- match.data(ps_matching)

library(survival)
library(gtsummary)

#Cox Regression for marginal HR
coxph(
  formula = Surv(eventtime, status) ~ exposure, 
  data = data_matched, 
  robust = TRUE, 
  weights = weights, 
  cluster = subclass
  ) |> 
  tbl_regression(exponentiate = TRUE)
```

<div id="dosfvpucgn" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#dosfvpucgn table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#dosfvpucgn thead, #dosfvpucgn tbody, #dosfvpucgn tfoot, #dosfvpucgn tr, #dosfvpucgn td, #dosfvpucgn th {
  border-style: none;
}
&#10;#dosfvpucgn p {
  margin: 0;
  padding: 0;
}
&#10;#dosfvpucgn .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#dosfvpucgn .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#dosfvpucgn .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#dosfvpucgn .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#dosfvpucgn .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#dosfvpucgn .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#dosfvpucgn .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#dosfvpucgn .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#dosfvpucgn .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#dosfvpucgn .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#dosfvpucgn .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#dosfvpucgn .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#dosfvpucgn .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#dosfvpucgn .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#dosfvpucgn .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#dosfvpucgn .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#dosfvpucgn .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#dosfvpucgn .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#dosfvpucgn .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#dosfvpucgn .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#dosfvpucgn .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#dosfvpucgn .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#dosfvpucgn .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#dosfvpucgn .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#dosfvpucgn .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#dosfvpucgn .gt_left {
  text-align: left;
}
&#10;#dosfvpucgn .gt_center {
  text-align: center;
}
&#10;#dosfvpucgn .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#dosfvpucgn .gt_font_normal {
  font-weight: normal;
}
&#10;#dosfvpucgn .gt_font_bold {
  font-weight: bold;
}
&#10;#dosfvpucgn .gt_font_italic {
  font-style: italic;
}
&#10;#dosfvpucgn .gt_super {
  font-size: 65%;
}
&#10;#dosfvpucgn .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#dosfvpucgn .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#dosfvpucgn .gt_indent_1 {
  text-indent: 5px;
}
&#10;#dosfvpucgn .gt_indent_2 {
  text-indent: 10px;
}
&#10;#dosfvpucgn .gt_indent_3 {
  text-indent: 15px;
}
&#10;#dosfvpucgn .gt_indent_4 {
  text-indent: 20px;
}
&#10;#dosfvpucgn .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;Characteristic&lt;/strong&gt;"><strong>Characteristic</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;HR&lt;/strong&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><strong>HR</strong><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;95% CI&lt;/strong&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><strong>95% CI</strong><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;p-value&lt;/strong&gt;"><strong>p-value</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">exposure</td>
<td headers="estimate" class="gt_row gt_center">0.82</td>
<td headers="ci" class="gt_row gt_center">0.75, 0.91</td>
<td headers="p.value" class="gt_row gt_center"><0.001</td></tr>
  </tbody>
  &#10;  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="4"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;"><sup>1</sup></span> HR = Hazard Ratio, CI = Confidence Interval</td>
    </tr>
  </tfoot>
</table>
</div>

## Appendix

### Repository details

This pre-populated example repository is conceptualized for a typical
RWE study with the following sub-directories/-files:

- `README.md`: Short description and information about the repository as
  well as mock study.

- `protocol`: HARPER quarto template (`protocol.qmd`) with empty table
  shells included in the protocol, template for study design figures and
  HARPER pdf template. It also contains a pre-configured
  `references.bib` BibTex file for citations/references one may want to
  include in the protocol.

- `documentation`: Sub-directory for documentation materials, e.g. data
  dictionaries or IRB approval letters.

- `scripts`: Sub-directory for R/Python analysis scripts. It is
  recommended to have a clear labeling such as: *1_data_query.qmd*,
  *2_descriptives.qmd*, *3_main_analysis.qmd*,
  *4_senstivity_analysis.qmd*, etc.

- `functions`: Custom R/Python functions to be sourced in analysis
  scripts.

- `results`: Sub-directory for collecting publication-ready tables,
  figures and other material relevant to communicate results.

- `public`: output directory for all quarto scripts and resources
  published on the website. This is particularly useful to publish
  annotated and interactive outputs along with a peer-reviewed
  manuscript as a measure of enhanced transparency and reproducibility
  and a way to augment and support study results. The name of the output
  directory can be customized in the `_quarto.yml` file according to
  specific system needs (for more details see the [quarto
  website](https://quarto.org/docs/projects/quarto-projects.html)).

- `manuscript`: Sub-directory where manuscript and supplementary
  material can be drafted and stored.

- `renv`: Project-specific R package library. Fore more information,
  please visit the [renv
  website](https://rstudio.github.io/renv/articles/renv.html).

- `.gitignore`: File to specify which files should not be tracked via
  git.

- `.Rprofile`: Environment file for R projects to store things like
  paths or keys (you can also setup an `.Renviron`)

- `.gitlab-ci.yml`: Pre-configured template. The `.gitlab-ci.yml`
  defines the GitLab continuous integration and deployment (CI/CD). This
  file may need to be adjusted to your specific project and installed
  GitLab runner(s).

Of course you can customize and remove/add other files and directories.
For example, one may also want to add bash files to orchestrate the
execution of different scripts.

## References

<div id="refs" class="references csl-bib-body">

<div id="ref-wang2022harmonized" class="csl-entry">

<span class="csl-left-margin">1.
</span><span class="csl-right-inline">Wang SV, Pottegård A, Crown W, *et
al.* HARmonized protocol template to enhance reproducibility of
hypothesis evaluating real-world evidence studies on treatment effects:
A good practices report of a joint ISPE/ISPOR task force. *Value in
Health* 2022; **25**: 1663–1672.</span>

</div>

<div id="ref-smdi" class="csl-entry">

<span class="csl-left-margin">2.
</span><span class="csl-right-inline">Weberpals J. Smdi: Perform
structural missing data investigations. 2023. Available at:
<https://janickweberpals.gitlab-pages.partners.org/smdi>.</span>

</div>

<div id="ref-MatchIt" class="csl-entry">

<span class="csl-left-margin">3.
</span><span class="csl-right-inline">Ho DE, Imai K, King G, Stuart EA.
MatchIt: Nonparametric preprocessing for parametric causal inference.
2011; **42**.
doi:[10.18637/jss.v042.i08](https://doi.org/10.18637/jss.v042.i08).</span>

</div>

</div>
