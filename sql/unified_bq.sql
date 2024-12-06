-- 1. Create the Combined DP Features Table
CREATE TABLE `rich-amphora-431916-m9.metadata_vcf_all.combined_dp_features_241029` AS

WITH first_features AS (
    SELECT
      RUN,
      COUNT(DP) AS total_variants,
      COUNTIF(G_AD_2 / DP > 0.5) AS major_variants,
      COUNTIF(G_AD_2 / DP < 0.5) AS minor_variants,
    
      -- Minor variant proportion of all variants
      COUNTIF(G_AD_2 / DP < 0.5) / COUNT(*) AS proportion_minor_variants,
    
      AVG(DP) AS mean_DP,
      STDDEV(DP) AS sd_DP,
    
      -- Raw fractions for different DP thresholds
      COUNTIF(DP < 10) / COUNT(*) AS fraction_below_10,
      COUNTIF(DP < 20) / COUNT(*) AS fraction_below_20,
      COUNTIF(DP < 30) / COUNT(*) AS fraction_below_30,
      COUNTIF(DP < 50) / COUNT(*) AS fraction_below_50,
      COUNTIF(DP < 100) / COUNT(*) AS fraction_below_100,
      COUNTIF(DP < 200) / COUNT(*) AS fraction_below_200,
      COUNTIF(DP < 300) / COUNT(*) AS fraction_below_300,
      COUNTIF(DP < 500) / COUNT(*) AS fraction_below_500,
      COUNTIF(DP < 1000) / COUNT(*) AS fraction_below_1000,
      COUNTIF(DP < 1500) / COUNT(*) AS fraction_below_1500,
      COUNTIF(DP < 2000) / COUNT(*) AS fraction_below_2000,
      COUNTIF(DP < 3000) / COUNT(*) AS fraction_below_3000,
      COUNTIF(DP < 5000) / COUNT(*) AS fraction_below_5000,
      COUNTIF(DP < 10000) / COUNT(*) AS fraction_below_10000,
    
      -- Log2 transformations for fractions with small offset
      LOG(COUNTIF(DP < 10) / COUNT(*) + 1e-4) / LOG(2) AS cs10,
      LOG(COUNTIF(DP < 20) / COUNT(*) + 1e-4) / LOG(2) AS cs20,
      LOG(COUNTIF(DP < 30) / COUNT(*) + 1e-4) / LOG(2) AS cs30,
      LOG(COUNTIF(DP < 50) / COUNT(*) + 1e-4) / LOG(2) AS cs50,
      LOG(COUNTIF(DP < 100) / COUNT(*) + 1e-4) / LOG(2) AS cs100,
      LOG(COUNTIF(DP < 200) / COUNT(*) + 1e-4) / LOG(2) AS cs200,
      LOG(COUNTIF(DP < 300) / COUNT(*) + 1e-4) / LOG(2) AS cs300,
      LOG(COUNTIF(DP < 500) / COUNT(*) + 1e-4) / LOG(2) AS cs500,
      LOG(COUNTIF(DP < 1000) / COUNT(*) + 1e-4) / LOG(2) AS cs1000,
      LOG(COUNTIF(DP < 1500) / COUNT(*) + 1e-4) / LOG(2) AS cs1500,
      LOG(COUNTIF(DP < 2000) / COUNT(*) + 1e-4) / LOG(2) AS cs2000,
      LOG(COUNTIF(DP < 3000) / COUNT(*) + 1e-4) / LOG(2) AS cs3000,
      LOG(COUNTIF(DP < 5000) / COUNT(*) + 1e-4) / LOG(2) AS cs5000,
      LOG(COUNTIF(DP < 10000) / COUNT(*) + 1e-4) / LOG(2) AS cs10000,
    
      -- Chaos calculation: -mean_DP / sd_DP
      -AVG(DP) / STDDEV(DP) AS chaos
    FROM
      `bigquery-public-data.nih_sequence_read.annotated_variations`
    GROUP BY
      RUN
    HAVING
      COUNT(DP) >= 3  -- Ensure there are at least 3 DP values
      AND STDDEV(DP) != 0  -- Avoid division by zero for chaos calculation
),

second_features AS (
    WITH spans AS (
      SELECT 
        acc, 
        l_pos, 
        r_pos, 
        1 + r_pos - l_pos AS span, 
        l_bnd AS min_depth,
        LEAST(h_bnd, 10000) AS max_depth,  -- Cap h_bnd at 10,000 to avoid excessively large values
        
        -- Compute bin width using capped max_depth
        GREATEST(LEAST(h_bnd, 10000) - l_bnd, 1) AS bin_width,
        
        -- Compute bin_weight using l_bnd, exclude infinite bin
        CASE
          WHEN h_bnd = 99999999 THEN NULL
          ELSE l_bnd
        END AS bin_weight,
        
        -- Compute log-weighted bin width (log scale) using l_bnd, exclude infinite bin
        CASE
          WHEN h_bnd = 99999999 THEN NULL
          ELSE 1.0 / LOG10(GREATEST(l_bnd, 1) + 1)  -- Adding 1 to avoid log(0) and handle l_bnd=0
        END AS bin_weight_log
      FROM 
        `bigquery-public-data.nih_sequence_read.sc2_depth_bins`
    ),
    stats AS (
      SELECT
        acc AS RUN,
        
        -- Count of spans per RUN
        COUNT(*) AS span_count,
        
        -- Weighted span count using l_bnd (exclude infinite bin)
        SUM(bin_weight) AS weighted_span_count,
        
        -- Log-weighted span count using l_bnd (exclude infinite bin)
        SUM(bin_weight_log) AS log_weighted_span_count,
        
        -- Total genome length
        MAX(r_pos) AS total_length,
        
        -- Weighted sums for mean and variance calculations using min_depth
        SUM(span * min_depth) AS sum_span_min_depth,
        SUM(span * POW(min_depth, 2)) AS sum_span_min_depth_sq,
    
        -- Cumulative fractions up to each bin's upper boundary
        SUM(CASE WHEN max_depth <= 0 THEN span ELSE 0 END) / MAX(r_pos) AS fraction_below_1,
        SUM(CASE WHEN max_depth <= 9 THEN span ELSE 0 END) / MAX(r_pos) AS fraction_below_10,
        SUM(CASE WHEN max_depth <= 49 THEN span ELSE 0 END) / MAX(r_pos) AS fraction_below_50,
        SUM(CASE WHEN max_depth <= 99 THEN span ELSE 0 END) / MAX(r_pos) AS fraction_below_100,
        SUM(CASE WHEN max_depth <= 999 THEN span ELSE 0 END) / MAX(r_pos) AS fraction_below_1000
      FROM
        spans
      GROUP BY
        RUN
    ),
    final_stats AS (
      SELECT
        RUN,
        span_count,
        weighted_span_count,
        log_weighted_span_count,
        total_length,
        sum_span_min_depth,
        sum_span_min_depth_sq,
        
        -- Mean depth calculation
        sum_span_min_depth / total_length AS mean_DP,
        
        -- Variance calculation using the weighted formula
        (sum_span_min_depth_sq / total_length) - POW(sum_span_min_depth / total_length, 2) AS variance_DP,
        
        -- Fractions at bin boundaries
        fraction_below_1,
        fraction_below_10,
        fraction_below_50,
        fraction_below_100,
        fraction_below_1000
      FROM
        stats
    ),
    computed_stats AS (
      SELECT
        RUN,
        span_count AS span_span_count,               -- Original span count
        weighted_span_count AS span_weighted_span_count,      -- Weighted span count using l_bnd
        log_weighted_span_count AS span_log_weighted_span_count,  -- Log-weighted span count using l_bnd
        mean_DP AS span_mean_DP,
        IF(variance_DP > 0, SQRT(variance_DP), NULL) AS span_sd_DP,
    
        -- Fractions at bin boundaries
        fraction_below_1 AS span_fraction_below_1,
        fraction_below_10 AS span_fraction_below_10,
        fraction_below_50 AS span_fraction_below_50,
        fraction_below_100 AS span_fraction_below_100,
        fraction_below_1000 AS span_fraction_below_1000,
    
        -- Log2 transformations for fractions with small offset
        LOG(fraction_below_1 + 1e-4) / LOG(2) AS span_cs1,
        LOG(fraction_below_10 + 1e-4) / LOG(2) AS span_cs10,
        LOG(fraction_below_50 + 1e-4) / LOG(2) AS span_cs50,
        LOG(fraction_below_100 + 1e-4) / LOG(2) AS span_cs100,
        LOG(fraction_below_1000 + 1e-4) / LOG(2) AS span_cs1000,
    
        -- Chaos calculation
        CASE
          WHEN variance_DP > 0 THEN -mean_DP / SQRT(variance_DP)
          ELSE NULL
        END AS span_chaos
      FROM
        final_stats
    )
    SELECT
      *
    FROM
      computed_stats
)

SELECT
    t1.RUN,
    
    -- Columns from the first_features
    t1.total_variants,
    t1.major_variants,
    t1.minor_variants,
    t1.proportion_minor_variants,
    t1.mean_DP AS dp_mean_DP,
    t1.sd_DP AS dp_sd_DP,
    t1.fraction_below_10,
    t1.fraction_below_20,
    t1.fraction_below_30,
    t1.fraction_below_50,
    t1.fraction_below_100,
    t1.fraction_below_200,
    t1.fraction_below_300,
    t1.fraction_below_500,
    t1.fraction_below_1000,
    t1.fraction_below_1500,
    t1.fraction_below_2000,
    t1.fraction_below_3000,
    t1.fraction_below_5000,
    t1.fraction_below_10000,
    t1.cs10,
    t1.cs20,
    t1.cs30,
    t1.cs50,
    t1.cs100,
    t1.cs200,
    t1.cs300,
    t1.cs500,
    t1.cs1000,
    t1.cs1500,
    t1.cs2000,
    t1.cs3000,
    t1.cs5000,
    t1.cs10000,
    t1.chaos AS dp_chaos,
    
    -- Columns from the second_features
    t2.span_span_count,
    t2.span_weighted_span_count,
    t2.span_log_weighted_span_count,
    t2.span_mean_DP,
    t2.span_sd_DP,
    t2.span_fraction_below_1,
    t2.span_fraction_below_10,
    t2.span_fraction_below_50,
    t2.span_fraction_below_100,
    t2.span_fraction_below_1000,
    t2.span_cs1,
    t2.span_cs10,
    t2.span_cs50,
    t2.span_cs100,
    t2.span_cs1000,
    t2.span_chaos
FROM
    first_features AS t1
INNER JOIN
    second_features AS t2
ON
    t1.RUN = t2.RUN;


-- 2. Create the Flat Metadata Table without Attributes
CREATE TABLE `rich-amphora-431916-m9.metadata_vcf_all.meta_flat_oct2024` AS 
SELECT DISTINCT
    metadata.* EXCEPT (
        jattr, 
        attributes, 
        datastore_filetype, 
        datastore_provider, 
        datastore_region,
        geo_loc_name_sam, 
        ena_first_public_run, 
        ena_last_update_run, 
        sample_name_sam, 
        biosamplemodel_sam
    )
FROM `nih-sra-datastore.sra.metadata` AS metadata
INNER JOIN `rich-amphora-431916-m9.metadata_vcf_all.combined_dp_features_241029` AS dp_features
    ON metadata.acc = dp_features.RUN
WHERE NOT EXISTS (
    SELECT 1
    FROM UNNEST(metadata.attributes) AS attribute
    WHERE REGEXP_CONTAINS(LOWER(attribute.k), r'(ww|sludge|wastewater|sewage)')
       OR REGEXP_CONTAINS(LOWER(attribute.v), r'(ww|sludge|wastewater|sewage)')
);

-- 3. Create the Attributes Only Table
CREATE TABLE `rich-amphora-431916-m9.metadata_vcf_all.attributes_only_oct2024` AS
SELECT 
    metadata.acc,
    attribute.k AS attribute_key,
    attribute.v AS attribute_value
FROM `nih-sra-datastore.sra.metadata` AS metadata
INNER JOIN `rich-amphora-431916-m9.metadata_vcf_all.combined_dp_features_241029` AS dp_features
    ON metadata.acc = dp_features.RUN
LEFT JOIN UNNEST(metadata.attributes) AS attribute
WHERE NOT EXISTS (
    SELECT 1
    FROM UNNEST(metadata.attributes) AS attr_filter
    WHERE REGEXP_CONTAINS(LOWER(attr_filter.k), r'(ww|sludge|wastewater|sewage)')
       OR REGEXP_CONTAINS(LOWER(attr_filter.v), r'(ww|sludge|wastewater|sewage)')
);

-- 4. Create the CT Only Table
CREATE TABLE `rich-amphora-431916-m9.metadata_vcf_all.ct_only_oct2024` AS
SELECT 
    acc AS sample_name, 
    attribute_key, 
    attribute_value
FROM `rich-amphora-431916-m9.metadata_vcf_all.attributes_only_oct2024` AS meta
WHERE 
    (
        REGEXP_CONTAINS(LOWER(meta.attribute_key), r'\bct\b')  -- 'ct' as a standalone word
        OR LOWER(meta.attribute_key) LIKE '%pcr%' 
        OR LOWER(meta.attribute_key) LIKE '%diag%'
    )
    AND REGEXP_CONTAINS(meta.attribute_value, r'^\d+(\.\d+)?$')  -- Checks if the value contains only numbers (with optional decimal)
    AND SAFE_CAST(meta.attribute_value AS FLOAT64) BETWEEN 5 AND 50
    AND LOWER(meta.attribute_value) NOT LIKE '%vol%'  -- Exclude if 'vol' is in the value
    AND LOWER(meta.attribute_value) NOT LIKE '%ww%';  -- Exclude if 'ww' is in the value

-- 5. Create the VCF Age Table
CREATE TABLE `rich-amphora-431916-m9.metadata_vcf_all.vcf_age_oct_2024` AS
SELECT 
    acc AS sample_name,
    attribute_key,
    attribute_value
FROM `rich-amphora-431916-m9.metadata_vcf_all.age_oct_2024` AS vax
INNER JOIN `rich-amphora-431916-m9.metadata_vcf_all.combined_dp_features_241029` AS dp
    ON vax.acc = dp.RUN
INNER JOIN `rich-amphora-431916-m9.metadata_vcf_all.meta_flat_oct2024` AS meta
    ON vax.acc = meta.acc;
