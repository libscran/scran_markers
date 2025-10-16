<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
  <compound kind="file">
    <name>score_markers_best.hpp</name>
    <path>/github/workspace/include/scran_markers/</path>
    <filename>score__markers__best_8hpp.html</filename>
    <class kind="struct">scran_markers::ScoreMarkersBestOptions</class>
    <class kind="struct">scran_markers::ScoreMarkersBestResults</class>
    <namespace>scran_markers</namespace>
  </compound>
  <compound kind="file">
    <name>score_markers_pairwise.hpp</name>
    <path>/github/workspace/include/scran_markers/</path>
    <filename>score__markers__pairwise_8hpp.html</filename>
    <class kind="struct">scran_markers::ScoreMarkersPairwiseOptions</class>
    <class kind="struct">scran_markers::ScoreMarkersPairwiseBuffers</class>
    <class kind="struct">scran_markers::ScoreMarkersPairwiseResults</class>
    <namespace>scran_markers</namespace>
  </compound>
  <compound kind="file">
    <name>score_markers_summary.hpp</name>
    <path>/github/workspace/include/scran_markers/</path>
    <filename>score__markers__summary_8hpp.html</filename>
    <includes id="summarize__comparisons_8hpp" name="summarize_comparisons.hpp" local="yes" import="no" module="no" objc="no">summarize_comparisons.hpp</includes>
    <class kind="struct">scran_markers::ScoreMarkersSummaryOptions</class>
    <class kind="struct">scran_markers::ScoreMarkersSummaryBuffers</class>
    <class kind="struct">scran_markers::ScoreMarkersSummaryResults</class>
    <namespace>scran_markers</namespace>
  </compound>
  <compound kind="file">
    <name>scran_markers.hpp</name>
    <path>/github/workspace/include/scran_markers/</path>
    <filename>scran__markers_8hpp.html</filename>
    <includes id="score__markers__pairwise_8hpp" name="score_markers_pairwise.hpp" local="yes" import="no" module="no" objc="no">score_markers_pairwise.hpp</includes>
    <includes id="score__markers__summary_8hpp" name="score_markers_summary.hpp" local="yes" import="no" module="no" objc="no">score_markers_summary.hpp</includes>
    <includes id="score__markers__best_8hpp" name="score_markers_best.hpp" local="yes" import="no" module="no" objc="no">score_markers_best.hpp</includes>
    <includes id="summarize__effects_8hpp" name="summarize_effects.hpp" local="yes" import="no" module="no" objc="no">summarize_effects.hpp</includes>
    <namespace>scran_markers</namespace>
  </compound>
  <compound kind="file">
    <name>summarize_comparisons.hpp</name>
    <path>/github/workspace/include/scran_markers/</path>
    <filename>summarize__comparisons_8hpp.html</filename>
    <class kind="struct">scran_markers::SummaryBuffers</class>
    <class kind="struct">scran_markers::SummaryResults</class>
    <namespace>scran_markers</namespace>
  </compound>
  <compound kind="file">
    <name>summarize_effects.hpp</name>
    <path>/github/workspace/include/scran_markers/</path>
    <filename>summarize__effects_8hpp.html</filename>
    <includes id="summarize__comparisons_8hpp" name="summarize_comparisons.hpp" local="yes" import="no" module="no" objc="no">summarize_comparisons.hpp</includes>
    <class kind="struct">scran_markers::SummarizeEffectsOptions</class>
    <namespace>scran_markers</namespace>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersBestOptions</name>
    <filename>structscran__markers_1_1ScoreMarkersBestOptions.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>threshold</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a1bfe79e22da0a7e1b8600baf13d3fbb3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>ad8fcee53b0fde5cf4e811746b96ef240</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_group_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a6c7a5f44323e9ce0e6910933a1cdcb30</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_group_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a69e2082f057cc0d2db9b66dbfc0635ea</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a123cba61d40fcdc36d8cc1e674b6965f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a2007c3ef0856d7c924820eeabe7c610f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a3c3663407f7f6634df71173b3e0415e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>aeb34e94fe7911ba89babdd338155d196</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>largest_cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a267157d5d176619cc7e0e6c322c1042d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>largest_auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>aa824707714a1690fbc23648578062e3b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>largest_delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a54def7db73a0f2eb9931443274b664d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>largest_delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>aa4a46e58c07d0c99483a254caeb3eba8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::optional&lt; double &gt;</type>
      <name>threshold_cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a929083a4a766df5e2703b40e135be433</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::optional&lt; double &gt;</type>
      <name>threshold_auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>acf10f836093794990810f8169c971ffa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::optional&lt; double &gt;</type>
      <name>threshold_delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a167e5d28ce2ea613a91fa96697834820</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::optional&lt; double &gt;</type>
      <name>threshold_delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>a94691cd88f06f8635561f85e1c70c0f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>keep_ties</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>aec6ceabe769286d92169913b28070ea9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::WeightPolicy</type>
      <name>block_weight_policy</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>ab5ef9efd63f9799d5d7023a8a2e8334c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::VariableWeightParameters</type>
      <name>variable_block_weight_parameters</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestOptions.html</anchorfile>
      <anchor>ae9155fca6deafc3acc0e4c92de471ea1</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersBestResults</name>
    <filename>structscran__markers_1_1ScoreMarkersBestResults.html</filename>
    <templarg>typename Stat_</templarg>
    <templarg>typename Index_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Stat_ &gt; &gt;</type>
      <name>mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestResults.html</anchorfile>
      <anchor>ad134d04290cc716dff38a32ba0282a6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Stat_ &gt; &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestResults.html</anchorfile>
      <anchor>af32e4f76f3572c9d686285c73e775bc8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; Index_, Stat_ &gt; &gt; &gt; &gt;</type>
      <name>cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestResults.html</anchorfile>
      <anchor>a7c6188aed4406851850793732aa030ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; Index_, Stat_ &gt; &gt; &gt; &gt;</type>
      <name>auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestResults.html</anchorfile>
      <anchor>a336708567f4dfacef8ba8bf188d720f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; Index_, Stat_ &gt; &gt; &gt; &gt;</type>
      <name>delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestResults.html</anchorfile>
      <anchor>a5cb20a7523c070686818ce889ea4d8ef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; Index_, Stat_ &gt; &gt; &gt; &gt;</type>
      <name>delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersBestResults.html</anchorfile>
      <anchor>a8a0be5bb72f393dddc153e40bb2cb8e2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersPairwiseBuffers</name>
    <filename>structscran__markers_1_1ScoreMarkersPairwiseBuffers.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Stat_ * &gt;</type>
      <name>mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseBuffers.html</anchorfile>
      <anchor>aff5d7b0d7873d3d4aa8dc00a03da3158</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ * &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseBuffers.html</anchorfile>
      <anchor>a44482d035e362711c472affb76289ea0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseBuffers.html</anchorfile>
      <anchor>a09a285a84e5ba22649d5aa0d87b2a21b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseBuffers.html</anchorfile>
      <anchor>ab052f7be5c0597a28f2372614a5bb3d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseBuffers.html</anchorfile>
      <anchor>a04fba49f4046f6ecc61edbc73ea3d546</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseBuffers.html</anchorfile>
      <anchor>a828364f6541a2c0ca63e9d05e020affc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersPairwiseOptions</name>
    <filename>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>threshold</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>a507dbb049c6a0412fbecebf245c77a57</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>a45f01e1456726ebf4211c969b5b4f128</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_group_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>a201e69fcbaacf2886b7056778f39a8a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_group_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>a933d78cc72b7cd25c91bc0aea014d1ad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>aafa7f898adc49ec6fdb72cc9300406c1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>ad4beffa1ddd3afabd66a3f086141b37e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>a4bcb5b1e56c1b50929961f59bd112577</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>af28ea3ee9dced2bfb608e23855a0fa0b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::WeightPolicy</type>
      <name>block_weight_policy</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>a24f56fd20bed4003ea896ee379735343</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::VariableWeightParameters</type>
      <name>variable_block_weight_parameters</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseOptions.html</anchorfile>
      <anchor>af6ab275bc541e6c0ffe62d9da271f550</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersPairwiseResults</name>
    <filename>structscran__markers_1_1ScoreMarkersPairwiseResults.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Stat_ &gt; &gt;</type>
      <name>mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseResults.html</anchorfile>
      <anchor>afd37e0bcf91b291101b96b40cad85675</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Stat_ &gt; &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseResults.html</anchorfile>
      <anchor>a32c8e8480dd9b4b67798613db3020632</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseResults.html</anchorfile>
      <anchor>a767e41fc513f9d35904267fe5376705b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseResults.html</anchorfile>
      <anchor>aeea05bc2c481d973297a7dfdb4409aeb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseResults.html</anchorfile>
      <anchor>a0bb5e0a19f2201c17069821234df873d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersPairwiseResults.html</anchorfile>
      <anchor>af80b58ddc75855f2a4f6117ed3b73cc0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersSummaryBuffers</name>
    <filename>structscran__markers_1_1ScoreMarkersSummaryBuffers.html</filename>
    <templarg>typename Stat_</templarg>
    <templarg>typename Rank_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Stat_ * &gt;</type>
      <name>mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryBuffers.html</anchorfile>
      <anchor>a84e6cc5b1e0ab7cff8e7c240f1b126cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ * &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryBuffers.html</anchorfile>
      <anchor>a77ee33709190fa5afced37804e62ae7f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryBuffers&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryBuffers.html</anchorfile>
      <anchor>ab2e2b1bde6d2b563ba95b91264a4538d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryBuffers&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryBuffers.html</anchorfile>
      <anchor>ad52aeadb45b041450c92ca4036ca39f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryBuffers&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryBuffers.html</anchorfile>
      <anchor>a9e2ed3d8b7a36bc5c8269a690f1e36c2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryBuffers&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryBuffers.html</anchorfile>
      <anchor>a7e5396d3a6530d20258ce45181e81ef9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersSummaryOptions</name>
    <filename>structscran__markers_1_1ScoreMarkersSummaryOptions.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>threshold</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a7b6daaa85845aac1cfa43b635ac43e9f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>acd9fd70b56bf779527de7a6068423bf6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_group_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a4ebd443ade8dc5f6f5d6602080201a59</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_group_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a1275da27456adda933f2516835199715</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>ab2100fe27f4471bdc296a718748ba699</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a4d614cc8e81337b5c5880a7adecdba7f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>ad10d15521459563cf74df5d7f9a98cc6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a2149d5a26bf1924de644cc0948f28839</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_min</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a51472efaed737fea2dc57f16fdc53870</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a76ada94f5fef79f959b41fd59c2d3224</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_median</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a660ccf44e64c70b72817ac3031b5c26b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_max</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a20a322581c10b9963c6ca639cd898cb0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_min_rank</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>abe5423cc6df71b83a0d6633979357f26</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::size_t</type>
      <name>min_rank_limit</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>af717200fc00125d5427c6f050fd50c89</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>min_rank_preserve_ties</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a865440ea066ee6a6cd81bfbbd9a8f426</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::WeightPolicy</type>
      <name>block_weight_policy</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a9581052b53623da99d9b1e6498d55430</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::VariableWeightParameters</type>
      <name>variable_block_weight_parameters</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>ad16d7561289a49c342ed6eb4be6adfbc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::ScoreMarkersSummaryResults</name>
    <filename>structscran__markers_1_1ScoreMarkersSummaryResults.html</filename>
    <templarg>typename Stat_</templarg>
    <templarg>typename Rank_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Stat_ &gt; &gt;</type>
      <name>mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryResults.html</anchorfile>
      <anchor>a1e5b0e717b37a567b014d5e9120704c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Stat_ &gt; &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryResults.html</anchorfile>
      <anchor>a80339ae0f037f4410bf03d4aed8cfefe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryResults&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>cohens_d</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryResults.html</anchorfile>
      <anchor>a19851939b8385ab173517b3e2e2acd49</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryResults&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>auc</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryResults.html</anchorfile>
      <anchor>a31be09bbdf2bf14b358895d1e455c1cf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryResults&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>delta_mean</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryResults.html</anchorfile>
      <anchor>aa5998e7e9bcd1a8c72145f7a8b63fd8b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; SummaryResults&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>delta_detected</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryResults.html</anchorfile>
      <anchor>a78e4984489cab29b9fd1cfeba65417b1</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::SummarizeEffectsOptions</name>
    <filename>structscran__markers_1_1SummarizeEffectsOptions.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__markers_1_1SummarizeEffectsOptions.html</anchorfile>
      <anchor>a11126e4aa2d7f2525d4ae98e4a2cf6dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_min</name>
      <anchorfile>structscran__markers_1_1SummarizeEffectsOptions.html</anchorfile>
      <anchor>aa882f8fbc4e85c2a0dae845d0a915169</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_mean</name>
      <anchorfile>structscran__markers_1_1SummarizeEffectsOptions.html</anchorfile>
      <anchor>a629404542bd27c11233b83f74b9cb30c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_median</name>
      <anchorfile>structscran__markers_1_1SummarizeEffectsOptions.html</anchorfile>
      <anchor>a263938e45059acbbb45c2aed3f1fea31</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_max</name>
      <anchorfile>structscran__markers_1_1SummarizeEffectsOptions.html</anchorfile>
      <anchor>aa006b2f98ae9a0b6bbbbd98bf2a497e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_min_rank</name>
      <anchorfile>structscran__markers_1_1SummarizeEffectsOptions.html</anchorfile>
      <anchor>a72ed641c77c6265d8da4583fe8f5db42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>min_rank_preserve_ties</name>
      <anchorfile>structscran__markers_1_1SummarizeEffectsOptions.html</anchorfile>
      <anchor>a5ea001f5ab4d49b472b55a9a79329230</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::SummaryBuffers</name>
    <filename>structscran__markers_1_1SummaryBuffers.html</filename>
    <templarg>typename Stat_</templarg>
    <templarg>typename Rank_</templarg>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>min</name>
      <anchorfile>structscran__markers_1_1SummaryBuffers.html</anchorfile>
      <anchor>aaa6f6973eed92f05f823376b9ecc326c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>mean</name>
      <anchorfile>structscran__markers_1_1SummaryBuffers.html</anchorfile>
      <anchor>a2f8f46dd9ba356d2e3da795445c4ef27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>median</name>
      <anchorfile>structscran__markers_1_1SummaryBuffers.html</anchorfile>
      <anchor>a8b9461e76a8d534bcf2cd53057c90d43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>max</name>
      <anchorfile>structscran__markers_1_1SummaryBuffers.html</anchorfile>
      <anchor>af471dad037118be3eda771d5a6b85b3e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Rank_ *</type>
      <name>min_rank</name>
      <anchorfile>structscran__markers_1_1SummaryBuffers.html</anchorfile>
      <anchor>afdf1258e4675f184d9caa7f4b5b283b7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_markers::SummaryResults</name>
    <filename>structscran__markers_1_1SummaryResults.html</filename>
    <templarg>typename Stat_</templarg>
    <templarg>typename Rank_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>min</name>
      <anchorfile>structscran__markers_1_1SummaryResults.html</anchorfile>
      <anchor>ab6661dc9ff5feafdfd4410493fd1d4e6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>mean</name>
      <anchorfile>structscran__markers_1_1SummaryResults.html</anchorfile>
      <anchor>a0b241991ebc25e59dfdf966e512c200a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>median</name>
      <anchorfile>structscran__markers_1_1SummaryResults.html</anchorfile>
      <anchor>a553f2a6628c8e5ce7a063a94b2be39fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>max</name>
      <anchorfile>structscran__markers_1_1SummaryResults.html</anchorfile>
      <anchor>aca9447a330e31fb89a04c00843a65a76</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Rank_ &gt;</type>
      <name>min_rank</name>
      <anchorfile>structscran__markers_1_1SummaryResults.html</anchorfile>
      <anchor>ab7ba8d989f18a8cc939f5c9435876800</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran_markers</name>
    <filename>namespacescran__markers.html</filename>
    <class kind="struct">scran_markers::ScoreMarkersBestOptions</class>
    <class kind="struct">scran_markers::ScoreMarkersBestResults</class>
    <class kind="struct">scran_markers::ScoreMarkersPairwiseBuffers</class>
    <class kind="struct">scran_markers::ScoreMarkersPairwiseOptions</class>
    <class kind="struct">scran_markers::ScoreMarkersPairwiseResults</class>
    <class kind="struct">scran_markers::ScoreMarkersSummaryBuffers</class>
    <class kind="struct">scran_markers::ScoreMarkersSummaryOptions</class>
    <class kind="struct">scran_markers::ScoreMarkersSummaryResults</class>
    <class kind="struct">scran_markers::SummarizeEffectsOptions</class>
    <class kind="struct">scran_markers::SummaryBuffers</class>
    <class kind="struct">scran_markers::SummaryResults</class>
    <member kind="function">
      <type>void</type>
      <name>score_markers_pairwise</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a6c63a08774d6f2105b68b656f8a3da94</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const ScoreMarkersPairwiseOptions &amp;options, const ScoreMarkersPairwiseBuffers&lt; Stat_ &gt; &amp;output)</arglist>
      <docanchor file="namespacescran__markers.html" title="Choice of effect size">effect-sizes</docanchor>
      <docanchor file="namespacescran__markers.html" title="With a minimum change threshold">threshold</docanchor>
      <docanchor file="namespacescran__markers.html" title="Other statistics">other</docanchor>
    </member>
    <member kind="function">
      <type>void</type>
      <name>score_markers_pairwise_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a9c5a2571c6f34d150a9b155ebd396a25</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const Block_ *const block, const ScoreMarkersPairwiseOptions &amp;options, const ScoreMarkersPairwiseBuffers&lt; Stat_ &gt; &amp;output)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersPairwiseResults&lt; Stat_ &gt;</type>
      <name>score_markers_pairwise</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>aa0e6a1d4dc2f4bad9e80ebe68eb08ce1</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const ScoreMarkersPairwiseOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersPairwiseResults&lt; Stat_ &gt;</type>
      <name>score_markers_pairwise_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>ab07e111b7e98525a9230937bb5ff61f8</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const Block_ *const block, const ScoreMarkersPairwiseOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersBestResults&lt; Stat_, Index_ &gt;</type>
      <name>score_markers_best</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a3b877db1571dd74a128af35eee5641d3</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, int top, const ScoreMarkersBestOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersBestResults&lt; Stat_, Index_ &gt;</type>
      <name>score_markers_best_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>afb11bdcf5c54fc4886a20e5f2f003aa2</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const Block_ *const block, int top, const ScoreMarkersBestOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>summarize_effects</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>aefddf4357cc7de227c5242b3861d5a27</anchor>
      <arglist>(const Gene_ ngenes, const std::size_t ngroups, const Stat_ *const effects, const std::vector&lt; SummaryBuffers&lt; Stat_, Rank_ &gt; &gt; &amp;summaries, const SummarizeEffectsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; SummaryResults&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>summarize_effects</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a5ae49bdf2cb9840a8dddf919de36107f</anchor>
      <arglist>(const Gene_ ngenes, const std::size_t ngroups, const Stat_ *const effects, const SummarizeEffectsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>score_markers_summary</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a92c28687f963fbe5535e3fc0c4a60a80</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const ScoreMarkersSummaryOptions &amp;options, const ScoreMarkersSummaryBuffers&lt; Stat_, Rank_ &gt; &amp;output)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>score_markers_summary_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>ab817d3509bf11224c5841cee6eb19d26</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const Block_ *const block, const ScoreMarkersSummaryOptions &amp;options, const ScoreMarkersSummaryBuffers&lt; Stat_, Rank_ &gt; &amp;output)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersSummaryResults&lt; Stat_, Rank_ &gt;</type>
      <name>score_markers_summary</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a2c9b15716f2c1a1d9d9add3d12319c7b</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const ScoreMarkersSummaryOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersSummaryResults&lt; Stat_, Rank_ &gt;</type>
      <name>score_markers_summary_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a3c5621dc5ac02a137bee6fab5c1c84e0</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *const group, const Block_ *const block, const ScoreMarkersSummaryOptions &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Marker detection for groups of cells</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
