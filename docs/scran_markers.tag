<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
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
      <type>int</type>
      <name>cache_size</name>
      <anchorfile>structscran__markers_1_1ScoreMarkersSummaryOptions.html</anchorfile>
      <anchor>a29fbb41c392819f973c9a96826555aa9</anchor>
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
      <anchor>a948a5ee108d7c1d07b9fec36553426e3</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const ScoreMarkersPairwiseOptions &amp;options, const ScoreMarkersPairwiseBuffers&lt; Stat_ &gt; &amp;output)</arglist>
      <docanchor file="namespacescran__markers.html" title="Choice of effect sizes">effect-sizes</docanchor>
      <docanchor file="namespacescran__markers.html" title="With a minimum change threshold">threshold</docanchor>
      <docanchor file="namespacescran__markers.html" title="Other statistics">other</docanchor>
    </member>
    <member kind="function">
      <type>void</type>
      <name>score_markers_pairwise_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a996c8a196ffe2f391028c9d33636e746</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const Block_ *block, const ScoreMarkersPairwiseOptions &amp;options, const ScoreMarkersPairwiseBuffers&lt; Stat_ &gt; &amp;output)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersPairwiseResults&lt; Stat_ &gt;</type>
      <name>score_markers_pairwise</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a2cbb54ae39a763ed86efa2f0a98f4482</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const ScoreMarkersPairwiseOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersPairwiseResults&lt; Stat_ &gt;</type>
      <name>score_markers_pairwise_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a08af688c64648170bc4e821a116c3eb8</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const Block_ *block, const ScoreMarkersPairwiseOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>summarize_effects</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a2a695132f3b280be336a1468b234f8b2</anchor>
      <arglist>(Gene_ ngenes, std::size_t ngroups, const Stat_ *effects, const std::vector&lt; SummaryBuffers&lt; Stat_, Rank_ &gt; &gt; &amp;summaries, const SummarizeEffectsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; SummaryResults&lt; Stat_, Rank_ &gt; &gt;</type>
      <name>summarize_effects</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a23e676487add3fc22adb886f42a49260</anchor>
      <arglist>(Gene_ ngenes, std::size_t ngroups, const Stat_ *effects, const SummarizeEffectsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>score_markers_summary</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a096a3c240fe439c2a47b9d83be2e50f6</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const ScoreMarkersSummaryOptions &amp;options, const ScoreMarkersSummaryBuffers&lt; Stat_, Rank_ &gt; &amp;output)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>score_markers_summary_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>acf8ecb4e2a494e5b442fe74cb30d9d92</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const Block_ *block, const ScoreMarkersSummaryOptions &amp;options, const ScoreMarkersSummaryBuffers&lt; Stat_, Rank_ &gt; &amp;output)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersSummaryResults&lt; Stat_, Rank_ &gt;</type>
      <name>score_markers_summary</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>a71bd424c39433e66942d0a89ae3b405d</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const ScoreMarkersSummaryOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMarkersSummaryResults&lt; Stat_, Rank_ &gt;</type>
      <name>score_markers_summary_blocked</name>
      <anchorfile>namespacescran__markers.html</anchorfile>
      <anchor>abb8febc2bb7932101ea0a9d4ab3d51f4</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Group_ *group, const Block_ *block, const ScoreMarkersSummaryOptions &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Marker detection for groups of cells</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
