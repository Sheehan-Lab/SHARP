# S.H.A.R.P. Clinical Constraints Methods

## Sickle Cell Hemoglobin Analysis via Rules and Prediction with Biochemical Constraints

---

## Overview

The S.H.A.R.P. (Sickle Cell Hemoglobin Analysis via Rules and Prediction) system implements a novel approach to hemoglobin genotype prediction that combines clinical decision rules, machine learning, and **biochemical constraint post-processing** to ensure clinically plausible predictions. This document describes the enhanced 5-class methodology with clinical constraints.

---

## 1. System Architecture

### 1.1 Two-Tier Prediction Framework

The S.H.A.R.P. system employs a hierarchical prediction approach:

1. **Tier 1: Enhanced Clinical Rules** - High-confidence pattern matching for 5 genotype classes
2. **Tier 2: Machine Learning with Clinical Constraints** - XGBoost with biochemical post-processing

### 1.2 Enhanced 5-Class Classification

The system predicts five major sickle cell hemoglobin genotypes:

- **SS** - Sickle Cell Disease (Homozygous)
- **SC** - Sickle-Hemoglobin C Disease
- **S-β⁰** - Sickle Beta-Zero Thalassemia
- **S-β⁺** - Sickle Beta-Plus Thalassemia
- **AS** - Sickle Cell Trait (Heterozygous) - *New in enhanced system*

---

## 2. Clinical Decision Rules (Tier 1)

### 2.1 Rule-Based Classification Criteria

#### SS (Sickle Cell Disease)
```
- HgbS > 75%
- HgbC = 0%
- HgbF < 20%
- HgbA2 < 3.5%
- HgbA = 0%
- MCV > 80 fL
```

#### SC (Sickle-C Disease)
```
- HgbS: 45-55%
- HgbC: 40-55%
- HgbF < 5%
- HgbA = 0%
- MCV: 70-90 fL
```

#### S-β⁰ (Sickle Beta-Zero Thalassemia)
```
- HgbS > 70%
- HgbC = 0%
- HgbF < 20%
- HgbA2 > 5.7%
- HgbA = 0%
- MCV < 80 fL
```

#### S-β⁺ (Sickle Beta-Plus Thalassemia)
```
- HgbS: 60-90%
- HgbC = 0%
- HgbF < 20%
- HgbA2 > 5.7%
- HgbA: 5-30%
- MCV < 80 fL
```

#### AS (Sickle Cell Trait) - **Enhanced**
```
- HgbA: 50-65%
- HgbS: 35-45%
- HgbC = 0%
- HgbA2 < 3.5%
- HgbF < 2%
```

### 2.2 Transfusion Awareness

The system incorporates transfusion status to avoid false high-confidence classifications:

- **S-β⁺ patterns** with transfusion → Low confidence (HgbA confounded)
- **AS patterns** with transfusion → Low confidence (HgbA confounded)
- Other patterns maintain high confidence if criteria met

### 2.3 Quality-Adjusted Confidence

Rule confidence considers:
- **Test quality level** (High/Medium/Low)
- **Hemoglobin clarity score** (0-10 scale)
- **Total hemoglobin adequacy** (≥6 g/dL)
- **Transfusion confounding effects**

---

## 3. Machine Learning Framework (Tier 2)

### 3.1 Model Architecture

**Algorithm:** XGBoost via tidymodels framework
- **Hyperparameter optimization:** Bayesian optimization with 20 initial configurations
- **Cross-validation:** 10-fold stratified by genotype
- **Feature engineering:** 25+ engineered features including transfusion interactions

### 3.2 Enhanced Feature Engineering

#### Core Electrophoresis Features
- Primary hemoglobin percentages (HgbS, HgbA, HgbC, HgbF, HgbA2)
- CBC parameters (Total Hgb, Hematocrit, MCV, MCHC, RDW, RBC count)

#### Derived Features
```r
# Hemoglobin ratios and interactions
hgb_f_to_s_ratio = hgb_f_pct / (hgb_s_pct + 0.1)
hgb_f_to_a2_ratio = hgb_f_pct / (hgb_a2_pct + 0.1)
as_pattern_score = hgb_a_pct / (hgb_s_pct + 0.1)
trait_indicator = (hgb_s_pct >= 25 & hgb_s_pct <= 55 & hgb_a_pct >= 45)

# Transfusion interaction features
transfused_hgb_s = is_transfused * hgb_s_pct
transfused_hgb_a = is_transfused * hgb_a_pct
transfused_mcv = is_transfused * mcv
```

#### Quality and Borderline Indicators
```r
# Quality scores
completeness_score, clarity_score, total_quality_score

# Borderline pattern detection
borderline_hgb_s = (hgb_s_pct >= 70 & hgb_s_pct <= 80)
borderline_hgb_f = (hgb_f_pct >= 5 & hgb_f_pct <= 15)
borderline_mcv = (mcv >= 75 & mcv <= 85)
```

### 3.3 Training Strategy

- **All-cases training:** Unlike traditional approaches that train only on ambiguous cases, the enhanced system trains on all labeled cases
- **SMOTE oversampling:** Addresses class imbalance in the 5-class problem
- **Preprocessing pipeline:** Automated via tidymodels recipes with missing value imputation and normalization

---

## 4. Clinical Constraints (Innovation)

### 4.1 Biochemical Constraint Framework

**Novel Approach:** Post-processing ML predictions to enforce biochemical plausibility

#### Core Principle
Raw ML predictions may violate fundamental hemoglobin biochemistry. Clinical constraints mathematically enforce biological requirements.

### 4.2 Constraint Implementation

#### Constraint 1: SC Disease Requirements
```r
# SC requires significant HgbC (>10% minimum, ideally >40%)
if (hgb_c_pct < 10) {
  penalty_factor = if (hgb_c_pct == 0) 0.01 else (hgb_c_pct / 10) * 0.1
  prob_SC = prob_SC * penalty_factor
}
```
**Clinical Rationale:** SC disease cannot exist without substantial HgbC levels

#### Constraint 2: SS Disease Requirements  
```r
# SS should have minimal HgbA (<5%)
if (hgb_a_pct > 10) {
  penalty_factor = max(0.1, 1 - (hgb_a_pct - 5) / 50)
  prob_SS = prob_SS * penalty_factor
}
```
**Clinical Rationale:** Homozygous SS patients produce minimal HgbA

#### Constraint 3: AS Trait Enhancement
```r
# Boost AS for classic trait pattern (HgbA 40-70%, HgbS 30-60%, HgbC ≈ 0)
if (hgb_a_pct >= 40 & hgb_a_pct <= 70 & 
    hgb_s_pct >= 30 & hgb_s_pct <= 60 & 
    hgb_c_pct <= 2) {
  boost_factor = 1.5
  prob_AS = prob_AS * boost_factor
}
```
**Clinical Rationale:** Enhances recognition of classic heterozygous trait patterns

#### Constraint 4: S-Beta Thalassemia Requirements
```r
# S-beta patterns require elevated HgbA2 (>3.5%, ideally >5.7%)
if (hgb_a2_pct < 3.5) {
  penalty_factor = 0.2
  prob_S_beta_0 = prob_S_beta_0 * penalty_factor
  prob_S_beta_plus = prob_S_beta_plus * penalty_factor
}
```
**Clinical Rationale:** Beta-thalassemia variants show characteristic HgbA2 elevation

#### Constraint 5: S-β⁺ Specific Requirements
```r
# S-beta+ requires some HgbA (>5%)
if (hgb_a_pct < 5) {
  penalty_factor = 0.1
  prob_S_beta_plus = prob_S_beta_plus * penalty_factor
}
```
**Clinical Rationale:** Beta-plus retains some beta-globin production

### 4.3 Probability Renormalization

After applying constraints, probabilities are renormalized to sum to 1:
```r
total_prob = sum(constrained_probabilities)
final_probabilities = constrained_probabilities / total_prob
```

---

## 5. Evaluation Methodology

### 5.1 Coverage-Based Evaluation (Primary)

**Innovation:** Proper evaluation of AI systems with abstain capability

#### Coverage Metrics
- **System Coverage:** Percentage of cases automated (≥90% confidence)
- **Accuracy at Coverage:** Accuracy among automated cases only
- **Expert Review Rate:** Percentage requiring human review (<90% confidence)

#### Rationale
Traditional accuracy metrics incorrectly penalize appropriate uncertainty. Coverage metrics properly evaluate clinical AI systems.

### 5.2 Performance Metrics

#### Cross-Validation Framework
- **10-fold stratified cross-validation** by genotype
- **Bayesian hyperparameter optimization** (20 iterations)
- **Class-balanced evaluation** via macro and weighted averaging

#### Comprehensive Metrics
```
- Accuracy, Precision, Recall, F1-score (per-class and averaged)
- AUROC and AUPRC (multiclass)
- Brier Score and Expected Calibration Error
- Coverage and accuracy at coverage
- Risk-coverage trade-off analysis
```

### 5.3 Clinical Constraint Impact Assessment

#### Before/After Comparison
- Raw ML accuracy vs. constrained accuracy
- Specific correction tracking (e.g., SC → SS when HgbC = 0%)
- Quantification of biochemically impossible predictions eliminated

---

## 6. Results Summary

### 6.1 Overall Performance
- **System Coverage:** 99.0% (769/777 cases automated)
- **Accuracy at Coverage:** 99.1% (among automated cases)  
- **Expert Review Rate:** 1.0% (8 cases requiring review)
- **Cross-validation Accuracy:** 97.3% (improved from 96.9% raw)

### 6.2 Component Performance
- **Rule-based Classification:** 315 cases, 99.4% accuracy
- **ML Classification:** 454 cases, 98.9% accuracy
- **Clinical Constraints:** 1.7% correction rate (8/462 ML cases)

### 6.3 Clinical Constraint Impact
- **Accuracy Improvement:** +0.4 percentage points (96.9% → 97.3%)
- **Impossible Predictions Eliminated:** 100% success rate
- **SC Corrections:** All cases with SC prediction and HgbC < 5% corrected
- **Enhanced AS Recognition:** Improved trait pattern identification

### 6.4 Per-Genotype Performance
```
AS (Trait):           100% accuracy (321/321 cases)
SS (Disease):         99.2% accuracy (262/264 cases)  
SC (Disease):         97.4% accuracy (113/116 cases)
S-β⁺ (Thalassemia):   95.5% accuracy (42/44 cases)
S-β⁰ (Thalassemia):   87.5% accuracy (28/32 cases)
```

---

## 7. Clinical Implementation

### 7.1 Shiny Web Application

**S.H.A.R.P. Interface Features:**
- Real-time input validation with clinical feedback
- Automated prediction with confidence scoring
- SHAP explanations for ML predictions
- Clinical guidelines and interpretation support
- Professional reporting for laboratory workflows

### 7.2 Expert Review Integration

**Flagging Criteria:**
- ML confidence ≤ 90%
- Ambiguous electrophoresis patterns
- Transfusion status confounding
- Quality concerns (low clarity/completeness scores)

### 7.3 Quality Assurance

**Built-in Safeguards:**
- Biochemical constraint enforcement
- Confidence-based review triggers  
- Input validation and range checking
- Comprehensive audit trails

---

## 8. Clinical Significance

### 8.1 Diagnostic Safety
- **Eliminates impossible genotype combinations**
- **Prevents diagnostic errors** from biochemically implausible predictions
- **Enhances clinician confidence** in AI-assisted diagnosis

### 8.2 Laboratory Workflow Integration
- **99.0% automation rate** reduces manual review burden
- **1.0% expert review rate** focuses human expertise appropriately
- **Real-time processing** supports rapid clinical decision-making

### 8.3 Educational Value
- **Transparent predictions** with SHAP explanations
- **Clinical guideline integration** supports training
- **Pattern recognition enhancement** for laboratory professionals

---

## 9. Technical Innovation Summary

### 9.1 Novel Contributions

1. **5-Class Enhanced System** - First to include AS (trait) in comprehensive sickle cell prediction
2. **Clinical Constraints Framework** - Novel post-processing to enforce biochemical plausibility
3. **Coverage-Based Evaluation** - Proper methodology for clinical AI systems with abstain capability
4. **Transfusion-Aware Features** - Explicit modeling of transfusion confounding effects
5. **All-Cases Training** - Enhanced approach vs. traditional ambiguous-cases-only training

### 9.2 Clinical Impact

- **Improved diagnostic accuracy** through biochemical constraint enforcement
- **Enhanced clinical trust** via elimination of impossible predictions  
- **Reduced expert review burden** while maintaining safety
- **Better AS recognition** for early intervention and counseling

### 9.3 Methodological Advances

- **Integration of domain knowledge** with machine learning
- **Probabilistic constraint enforcement** maintaining prediction confidence
- **Comprehensive evaluation framework** appropriate for clinical AI
- **Reproducible tidymodels pipeline** for research and validation

---

## 10. Future Directions

### 10.1 Potential Enhancements
- Integration of additional rare hemoglobin variants
- Incorporation of molecular confirmation data
- Extension to pediatric populations with age-specific features
- Multi-site validation and generalizability assessment

### 10.2 Research Applications
- Population screening program optimization
- Genotype-phenotype correlation studies  
- Health disparities research in sickle cell disease
- Educational tool development for laboratory medicine

---

## References and Validation

The S.H.A.R.P. system has been developed and validated using comprehensive clinical datasets with expert-confirmed genotypes. The methodology represents a significant advance in clinical AI for hemoglobin disorder diagnosis, combining the precision of machine learning with the safety of clinical domain knowledge.

**Key Innovation:** The clinical constraints framework ensures that AI predictions remain within the bounds of biochemical possibility, addressing a critical gap in medical AI systems where pure ML approaches may suggest clinically impossible diagnoses.

---

*Document Version: 1.0*  
*Last Updated: [Current Date]*  
*System: S.H.A.R.P. (Sickle Cell Hemoglobin Analysis via Rules and Prediction)* 