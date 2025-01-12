# Set up

Set up workspace, set options, and load required packages.

    knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

Load libraries.

    library(tidyverse)
    library(ggplot2)
    library(survival)
    library(readxl)
    library(ggsurvfit)
    library(gtsummary)
    library(cardx)
    library(cowplot)

# Read data

Read in data.

    data<-read_excel("data/survival/survival_resazurin.xlsx")

Turn into long format.

    data<-data%>%
      pivot_longer(names_to="time", values_to="status", cols=`0`:`24`)%>%
      mutate(time=as.numeric(time))

Create data frames with control and high temperature subsets.

    control<-data%>%filter(temperature=="18C")
    stress<-data%>%filter(temperature=="42C")

# Generate Kaplan Meier survival curves

## Control temperature 18C

    s1 <- survfit(Surv(time, status) ~ hardening, data = control)
    str(s1)

    ## List of 18
    ##  $ n        : int [1:5] 360 240 240 240 240
    ##  $ time     : num [1:30] 0 1 2 3 4 24 0 1 2 3 ...
    ##  $ n.risk   : num [1:30] 360 300 240 180 120 60 240 200 160 120 ...
    ##  $ n.event  : num [1:30] 0 0 0 0 4 14 0 0 0 0 ...
    ##  $ n.censor : num [1:30] 60 60 60 60 56 46 40 40 40 40 ...
    ##  $ surv     : num [1:30] 1 1 1 1 0.967 ...
    ##  $ std.err  : num [1:30] 0 0 0 0 0.017 ...
    ##  $ cumhaz   : num [1:30] 0 0 0 0 0.0333 ...
    ##  $ std.chaz : num [1:30] 0 0 0 0 0.0167 ...
    ##  $ strata   : Named int [1:5] 6 6 6 6 6
    ##   ..- attr(*, "names")= chr [1:5] "hardening=control" "hardening=fresh-water" "hardening=fresh-water-temperature" "hardening=immune" ...
    ##  $ type     : chr "right"
    ##  $ logse    : logi TRUE
    ##  $ conf.int : num 0.95
    ##  $ conf.type: chr "log"
    ##  $ lower    : num [1:30] 1 1 1 1 0.935 ...
    ##  $ upper    : num [1:30] 1 1 1 1 0.999 ...
    ##  $ t0       : num 0
    ##  $ call     : language survfit(formula = Surv(time, status) ~ hardening, data = control)
    ##  - attr(*, "class")= chr "survfit"

Plot the survival function

    survfit2(Surv(time, status) ~ hardening+temperature, data = control) %>% 
      ggsurvfit() +
      labs(
        x = "Hours",
        y = "Survival probability"
      )

![](survival_files/figure-markdown_strict/unnamed-chunk-7-1.png)

Estimate survival at 24 hours.

    summary(survfit(Surv(time, status) ~ hardening, data = control), times = 24)

    ## Call: survfit(formula = Surv(time, status) ~ hardening, data = control)
    ## 
    ##                 hardening=control 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      60.0000      18.0000       0.7411       0.0543       0.6420 
    ## upper 95% CI 
    ##       0.8555 
    ## 
    ##                 hardening=fresh-water 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      40.0000      12.0000       0.7159       0.0703       0.5906 
    ## upper 95% CI 
    ##       0.8679 
    ## 
    ##                 hardening=fresh-water-temperature 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##       24.000       40.000       16.000        0.617        0.076        0.485 
    ## upper 95% CI 
    ##        0.786 
    ## 
    ##                 hardening=immune 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      40.0000      13.0000       0.7362       0.0655       0.6184 
    ## upper 95% CI 
    ##       0.8765 
    ## 
    ##                 hardening=temperature 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      40.0000      11.0000       0.7250       0.0706       0.5990 
    ## upper 95% CI 
    ##       0.8775

Estimated survival is ~87% in control temperature.

Use a log rank model to determine statistical differences in curves.

    survdiff(Surv(time, status) ~ hardening, data = control)

    ## Call:
    ## survdiff(formula = Surv(time, status) ~ hardening, data = control)
    ## 
    ##                                     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## hardening=control                 360       18     19.1   0.06234   0.11188
    ## hardening=fresh-water             240       12     12.7   0.04156   0.06630
    ## hardening=fresh-water-temperature 240       16     12.7   0.84156   1.34252
    ## hardening=immune                  240       13     12.7   0.00584   0.00932
    ## hardening=temperature             240       11     12.7   0.23442   0.37396
    ## 
    ##  Chisq= 1.5  on 4 degrees of freedom, p= 0.8

    # Call:
    # survdiff(formula = Surv(time, status) ~ hardening, data = control)
    # 
    #                                     N Observed Expected (O-E)^2/E (O-E)^2/V
    # hardening=control                 360       18     19.1   0.06234   0.11188
    # hardening=fresh-water             240       12     12.7   0.04156   0.06630
    # hardening=fresh-water-temperature 240       16     12.7   0.84156   1.34252
    # hardening=immune                  240       13     12.7   0.00584   0.00932
    # hardening=temperature             240       11     12.7   0.23442   0.37396
    # 
    #  Chisq= 1.5  on 4 degrees of freedom, p= 0.8 

There is no difference in survival at 18C at p=0.8.

Analyze again with a Cox proportional hazards model.

    coxph(Surv(time, status) ~ hardening, data = control)

    ## Call:
    ## coxph(formula = Surv(time, status) ~ hardening, data = control)
    ## 
    ##                                      coef exp(coef) se(coef)      z     p
    ## hardeningfresh-water              0.02005   1.02025  0.37269  0.054 0.957
    ## hardeningfresh-water-temperature  0.36099   1.43475  0.34368  1.050 0.294
    ## hardeningimmune                   0.07687   1.07991  0.36398  0.211 0.833
    ## hardeningtemperature             -0.06777   0.93447  0.38271 -0.177 0.859
    ## 
    ## Likelihood ratio test=1.57  on 4 df, p=0.8147
    ## n= 1320, number of events= 70

    # coxph(formula = Surv(time, status) ~ hardening, data = control)
    # 
    #                                      coef exp(coef) se(coef)      z     p
    # hardeningfresh-water              0.02005   1.02025  0.37269  0.054 0.957
    # hardeningfresh-water-temperature  0.36099   1.43475  0.34368  1.050 0.294
    # hardeningimmune                   0.07687   1.07991  0.36398  0.211 0.833
    # hardeningtemperature             -0.06777   0.93447  0.38271 -0.177 0.859
    # 
    # Likelihood ratio test=1.57  on 4 df, p=0.8147
    # n= 1320, number of events= 70  

    coxph(Surv(time, status) ~ hardening, data = control) %>% 
      tbl_regression(exp = TRUE) 

<div id="ncyyfblwqb" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ncyyfblwqb table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#ncyyfblwqb thead, #ncyyfblwqb tbody, #ncyyfblwqb tfoot, #ncyyfblwqb tr, #ncyyfblwqb td, #ncyyfblwqb th {
  border-style: none;
}

#ncyyfblwqb p {
  margin: 0;
  padding: 0;
}

#ncyyfblwqb .gt_table {
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

#ncyyfblwqb .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#ncyyfblwqb .gt_title {
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

#ncyyfblwqb .gt_subtitle {
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

#ncyyfblwqb .gt_heading {
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

#ncyyfblwqb .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ncyyfblwqb .gt_col_headings {
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

#ncyyfblwqb .gt_col_heading {
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

#ncyyfblwqb .gt_column_spanner_outer {
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

#ncyyfblwqb .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ncyyfblwqb .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ncyyfblwqb .gt_column_spanner {
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

#ncyyfblwqb .gt_spanner_row {
  border-bottom-style: hidden;
}

#ncyyfblwqb .gt_group_heading {
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

#ncyyfblwqb .gt_empty_group_heading {
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

#ncyyfblwqb .gt_from_md > :first-child {
  margin-top: 0;
}

#ncyyfblwqb .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ncyyfblwqb .gt_row {
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

#ncyyfblwqb .gt_stub {
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

#ncyyfblwqb .gt_stub_row_group {
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

#ncyyfblwqb .gt_row_group_first td {
  border-top-width: 2px;
}

#ncyyfblwqb .gt_row_group_first th {
  border-top-width: 2px;
}

#ncyyfblwqb .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ncyyfblwqb .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#ncyyfblwqb .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#ncyyfblwqb .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ncyyfblwqb .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ncyyfblwqb .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ncyyfblwqb .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#ncyyfblwqb .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ncyyfblwqb .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ncyyfblwqb .gt_footnotes {
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

#ncyyfblwqb .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ncyyfblwqb .gt_sourcenotes {
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

#ncyyfblwqb .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ncyyfblwqb .gt_left {
  text-align: left;
}

#ncyyfblwqb .gt_center {
  text-align: center;
}

#ncyyfblwqb .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ncyyfblwqb .gt_font_normal {
  font-weight: normal;
}

#ncyyfblwqb .gt_font_bold {
  font-weight: bold;
}

#ncyyfblwqb .gt_font_italic {
  font-style: italic;
}

#ncyyfblwqb .gt_super {
  font-size: 65%;
}

#ncyyfblwqb .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#ncyyfblwqb .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#ncyyfblwqb .gt_indent_1 {
  text-indent: 5px;
}

#ncyyfblwqb .gt_indent_2 {
  text-indent: 10px;
}

#ncyyfblwqb .gt_indent_3 {
  text-indent: 15px;
}

#ncyyfblwqb .gt_indent_4 {
  text-indent: 20px;
}

#ncyyfblwqb .gt_indent_5 {
  text-indent: 25px;
}

#ncyyfblwqb .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#ncyyfblwqb div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;Characteristic&lt;/strong&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>Characteristic</strong></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;HR&lt;/strong&gt;&lt;/span&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>HR</strong></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;95% CI&lt;/strong&gt;&lt;/span&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>95% CI</strong></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;p-value&lt;/strong&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>p-value</strong></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">hardening</td>
<td headers="estimate" class="gt_row gt_center"><br /></td>
<td headers="conf.low" class="gt_row gt_center"><br /></td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    control</td>
<td headers="estimate" class="gt_row gt_center">—</td>
<td headers="conf.low" class="gt_row gt_center">—</td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    fresh-water</td>
<td headers="estimate" class="gt_row gt_center">1.02</td>
<td headers="conf.low" class="gt_row gt_center">0.49, 2.12</td>
<td headers="p.value" class="gt_row gt_center">>0.9</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    fresh-water-temperature</td>
<td headers="estimate" class="gt_row gt_center">1.43</td>
<td headers="conf.low" class="gt_row gt_center">0.73, 2.81</td>
<td headers="p.value" class="gt_row gt_center">0.3</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    immune</td>
<td headers="estimate" class="gt_row gt_center">1.08</td>
<td headers="conf.low" class="gt_row gt_center">0.53, 2.20</td>
<td headers="p.value" class="gt_row gt_center">0.8</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    temperature</td>
<td headers="estimate" class="gt_row gt_center">0.93</td>
<td headers="conf.low" class="gt_row gt_center">0.44, 1.98</td>
<td headers="p.value" class="gt_row gt_center">0.9</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="4"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span> <span class='gt_from_md'>HR = Hazard Ratio, CI = Confidence Interval</span></td>
    </tr>
  </tfoot>
</table>
</div>

Control is the reference level. No differences in survival at 18C.

# Generate Kaplan Meier survival curves

## Stress temperature 42C

    s2 <- survfit(Surv(time, status) ~ hardening, data = stress)
    str(s2)

    ## List of 18
    ##  $ n        : int [1:5] 360 240 240 240 240
    ##  $ time     : num [1:30] 0 1 2 3 4 24 0 1 2 3 ...
    ##  $ n.risk   : num [1:30] 360 300 240 180 120 60 240 200 160 120 ...
    ##  $ n.event  : num [1:30] 0 1 14 21 32 48 0 1 12 16 ...
    ##  $ n.censor : num [1:30] 60 59 46 39 28 12 40 39 28 24 ...
    ##  $ surv     : num [1:30] 1 0.997 0.939 0.829 0.608 ...
    ##  $ std.err  : num [1:30] 0 0.00334 0.01641 0.03167 0.06351 ...
    ##  $ cumhaz   : num [1:30] 0 0.00333 0.06167 0.17833 0.445 ...
    ##  $ std.chaz : num [1:30] 0 0.00333 0.01594 0.03004 0.0559 ...
    ##  $ strata   : Named int [1:5] 6 6 6 6 6
    ##   ..- attr(*, "names")= chr [1:5] "hardening=control" "hardening=fresh-water" "hardening=fresh-water-temperature" "hardening=immune" ...
    ##  $ type     : chr "right"
    ##  $ logse    : logi TRUE
    ##  $ conf.int : num 0.95
    ##  $ conf.type: chr "log"
    ##  $ lower    : num [1:30] 1 0.99 0.909 0.779 0.537 ...
    ##  $ upper    : num [1:30] 1 1 0.969 0.882 0.689 ...
    ##  $ t0       : num 0
    ##  $ call     : language survfit(formula = Surv(time, status) ~ hardening, data = stress)
    ##  - attr(*, "class")= chr "survfit"

Plot the survival function

    survfit2(Surv(time, status) ~ hardening+temperature, data = stress) %>% 
      ggsurvfit() +
      labs(
        x = "Hours",
        y = "Survival probability"
      )

![](survival_files/figure-markdown_strict/unnamed-chunk-12-1.png)

Estimate survival at 24 hours.

    summary(survfit(Surv(time, status) ~ hardening, data = stress), times = 24)

    ## Call: survfit(formula = Surv(time, status) ~ hardening, data = stress)
    ## 
    ##                 hardening=control 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      60.0000     116.0000       0.1216       0.0323       0.0722 
    ## upper 95% CI 
    ##       0.2048 
    ## 
    ##                 hardening=fresh-water 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      40.0000      83.0000       0.1029       0.0363       0.0516 
    ## upper 95% CI 
    ##       0.2054 
    ## 
    ##                 hardening=fresh-water-temperature 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      40.0000      83.0000       0.1551       0.0418       0.0914 
    ## upper 95% CI 
    ##       0.2631 
    ## 
    ##                 hardening=immune 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      40.0000      97.0000       0.1386       0.0376       0.0814 
    ## upper 95% CI 
    ##       0.2360 
    ## 
    ##                 hardening=temperature 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##      24.0000      40.0000     109.0000       0.0549       0.0236       0.0236 
    ## upper 95% CI 
    ##       0.1277

Estimated stress survival is 12-26%. Temperature hardened is the lowest.

Use a log rank model to determine statistical differences in curves.

    survdiff(Surv(time, status) ~ hardening, data = stress)

    ## Call:
    ## survdiff(formula = Surv(time, status) ~ hardening, data = stress)
    ## 
    ##                                     N Observed Expected (O-E)^2/E (O-E)^2/V
    ## hardening=control                 360      116    133.1     2.195     5.048
    ## hardening=fresh-water             240       83     88.7     0.370     0.756
    ## hardening=fresh-water-temperature 240       83     88.7     0.370     0.756
    ## hardening=immune                  240       97     88.7     0.771     1.577
    ## hardening=temperature             240      109     88.7     4.632     9.470
    ## 
    ##  Chisq= 13.9  on 4 degrees of freedom, p= 0.007

    # Call:
    # survdiff(formula = Surv(time, status) ~ hardening, data = stress)
    # 
    #                                     N Observed Expected (O-E)^2/E (O-E)^2/V
    # hardening=control                 360      116    133.1     2.195     5.048
    # hardening=fresh-water             240       83     88.7     0.370     0.756
    # hardening=fresh-water-temperature 240       83     88.7     0.370     0.756
    # hardening=immune                  240       97     88.7     0.771     1.577
    # hardening=temperature             240      109     88.7     4.632     9.470
    # 
    #  Chisq= 13.9  on 4 degrees of freedom, p= 0.007 

There is a significant difference in survival at p=0.007.

Analyze again with a Cox proportional hazards model.

    coxph(Surv(time, status) ~ hardening, data = stress)

    ## Call:
    ## coxph(formula = Surv(time, status) ~ hardening, data = stress)
    ## 
    ##                                     coef exp(coef) se(coef)     z       p
    ## hardeningfresh-water             0.08431   1.08797  0.14377 0.586 0.55759
    ## hardeningfresh-water-temperature 0.04540   1.04644  0.14380 0.316 0.75224
    ## hardeningimmune                  0.20959   1.23317  0.13763 1.523 0.12780
    ## hardeningtemperature             0.40947   1.50602  0.13343 3.069 0.00215
    ## 
    ## Likelihood ratio test=11.11  on 4 df, p=0.02538
    ## n= 1320, number of events= 488

    # Call:
    # coxph(formula = Surv(time, status) ~ hardening, data = stress)
    # 
    #                                     coef exp(coef) se(coef)     z       p
    # hardeningfresh-water             0.08431   1.08797  0.14377 0.586 0.55759
    # hardeningfresh-water-temperature 0.04540   1.04644  0.14380 0.316 0.75224
    # hardeningimmune                  0.20959   1.23317  0.13763 1.523 0.12780
    # hardeningtemperature             0.40947   1.50602  0.13343 3.069 0.00215
    # 
    # Likelihood ratio test=11.11  on 4 df, p=0.02538
    # n= 1320, number of events= 488 

    coxph(Surv(time, status) ~ hardening, data = stress) %>% 
      tbl_regression(exp = TRUE) 

<div id="bprrdfvsar" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#bprrdfvsar table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#bprrdfvsar thead, #bprrdfvsar tbody, #bprrdfvsar tfoot, #bprrdfvsar tr, #bprrdfvsar td, #bprrdfvsar th {
  border-style: none;
}

#bprrdfvsar p {
  margin: 0;
  padding: 0;
}

#bprrdfvsar .gt_table {
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

#bprrdfvsar .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#bprrdfvsar .gt_title {
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

#bprrdfvsar .gt_subtitle {
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

#bprrdfvsar .gt_heading {
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

#bprrdfvsar .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bprrdfvsar .gt_col_headings {
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

#bprrdfvsar .gt_col_heading {
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

#bprrdfvsar .gt_column_spanner_outer {
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

#bprrdfvsar .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#bprrdfvsar .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#bprrdfvsar .gt_column_spanner {
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

#bprrdfvsar .gt_spanner_row {
  border-bottom-style: hidden;
}

#bprrdfvsar .gt_group_heading {
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

#bprrdfvsar .gt_empty_group_heading {
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

#bprrdfvsar .gt_from_md > :first-child {
  margin-top: 0;
}

#bprrdfvsar .gt_from_md > :last-child {
  margin-bottom: 0;
}

#bprrdfvsar .gt_row {
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

#bprrdfvsar .gt_stub {
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

#bprrdfvsar .gt_stub_row_group {
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

#bprrdfvsar .gt_row_group_first td {
  border-top-width: 2px;
}

#bprrdfvsar .gt_row_group_first th {
  border-top-width: 2px;
}

#bprrdfvsar .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bprrdfvsar .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#bprrdfvsar .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#bprrdfvsar .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bprrdfvsar .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bprrdfvsar .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#bprrdfvsar .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#bprrdfvsar .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#bprrdfvsar .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bprrdfvsar .gt_footnotes {
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

#bprrdfvsar .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bprrdfvsar .gt_sourcenotes {
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

#bprrdfvsar .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bprrdfvsar .gt_left {
  text-align: left;
}

#bprrdfvsar .gt_center {
  text-align: center;
}

#bprrdfvsar .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#bprrdfvsar .gt_font_normal {
  font-weight: normal;
}

#bprrdfvsar .gt_font_bold {
  font-weight: bold;
}

#bprrdfvsar .gt_font_italic {
  font-style: italic;
}

#bprrdfvsar .gt_super {
  font-size: 65%;
}

#bprrdfvsar .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#bprrdfvsar .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#bprrdfvsar .gt_indent_1 {
  text-indent: 5px;
}

#bprrdfvsar .gt_indent_2 {
  text-indent: 10px;
}

#bprrdfvsar .gt_indent_3 {
  text-indent: 15px;
}

#bprrdfvsar .gt_indent_4 {
  text-indent: 20px;
}

#bprrdfvsar .gt_indent_5 {
  text-indent: 25px;
}

#bprrdfvsar .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#bprrdfvsar div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;Characteristic&lt;/strong&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>Characteristic</strong></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;HR&lt;/strong&gt;&lt;/span&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>HR</strong></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;95% CI&lt;/strong&gt;&lt;/span&gt;&lt;span class=&quot;gt_footnote_marks&quot; style=&quot;white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;&quot;&gt;&lt;sup&gt;1&lt;/sup&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>95% CI</strong></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;span class='gt_from_md'&gt;&lt;strong&gt;p-value&lt;/strong&gt;&lt;/span&gt;"><span class='gt_from_md'><strong>p-value</strong></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">hardening</td>
<td headers="estimate" class="gt_row gt_center"><br /></td>
<td headers="conf.low" class="gt_row gt_center"><br /></td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    control</td>
<td headers="estimate" class="gt_row gt_center">—</td>
<td headers="conf.low" class="gt_row gt_center">—</td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    fresh-water</td>
<td headers="estimate" class="gt_row gt_center">1.09</td>
<td headers="conf.low" class="gt_row gt_center">0.82, 1.44</td>
<td headers="p.value" class="gt_row gt_center">0.6</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    fresh-water-temperature</td>
<td headers="estimate" class="gt_row gt_center">1.05</td>
<td headers="conf.low" class="gt_row gt_center">0.79, 1.39</td>
<td headers="p.value" class="gt_row gt_center">0.8</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    immune</td>
<td headers="estimate" class="gt_row gt_center">1.23</td>
<td headers="conf.low" class="gt_row gt_center">0.94, 1.62</td>
<td headers="p.value" class="gt_row gt_center">0.13</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    temperature</td>
<td headers="estimate" class="gt_row gt_center">1.51</td>
<td headers="conf.low" class="gt_row gt_center">1.16, 1.96</td>
<td headers="p.value" class="gt_row gt_center">0.002</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="4"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height: 0;"><sup>1</sup></span> <span class='gt_from_md'>HR = Hazard Ratio, CI = Confidence Interval</span></td>
    </tr>
  </tfoot>
</table>
</div>

Control is the reference level. There is a significant difference
between control and temperature hardening survival.

# Generate plots

Set theme.

    my_theme<-theme_classic()

Control

    plot1<-survfit2(Surv(time, status) ~ hardening, data = control) %>% 
      ggsurvfit() +
      labs(
        x = "Hours",
        y = "Survival probability",
        title="18°C"
      )+
      ylim(0,1)+
      my_theme+
      geom_text(x=10, y=0.2, label="Cox PH p=0.8")+
      theme(legend.position="none")

Stress

    plot2<-survfit2(Surv(time, status) ~ hardening, data = stress) %>% 
      ggsurvfit() +
      labs(
        x = "Hours",
        y = "Survival probability",
        title = "42°C"
      )+
      ylim(0,1)+
      my_theme+
      geom_text(x=10, y=0.2, label="Cox PH p=0.025")+
      geom_text(x=10, y=0.15, label="Temperature vs control p=0.002")+
      geom_text(x=10, y=0.10, label="All others vs control p>0.05")+
      theme(legend.position="right")

Assemble plot

    plots<-plot_grid(plot1, plot2, rel_widths=c(0.65,1), ncol=2)

    ggsave(plots, filename="figures/survival/KMcurves.png", width=10, height=4)
