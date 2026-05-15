# Configuration Options

UKB-GYM applets are highly customizable through YAML configuration files found in the `configs/` directory. While the applets come bundled with default configurations, you can override these defaults by providing your own YAML files via the CLI or GUI.

## 1. config_variant_classes.yaml

This file defines the preset variant filtering classes available across all applets. When you set the `-i variant_class="missense"` input parameter, the applet looks up the "missense" definition in this config.

**Structure:**
- Keys represent the `variant_class` name (e.g., `missense`, `intron`, `enhancer_encode`).
- `variant_filtering`: A list of string expressions (written in Polars filter syntax) used to dynamically subset the variant annotations dataframe.
- `tool_categories`: A list specifying which annotation categories are relevant to this class.
- `x_label`: The display name used in plot generation.

**Custom Override:**
If you need to analyze a completely novel class of variants, you can set `-i variant_class="other"` and provide a custom YAML file mapping your filters via the `variant_class_config_yaml` input.

## 2. config_correlations.yaml

Used exclusively by the `correlations_traits` and `correlations_olink` applets. This file defines the visual metadata for each scoring tool when plotting correlation results.

**Structure:**
- Keys represent annotation categories (e.g., `plof`, `missense`, `conservation`).
- Inside each category, individual tools are mapped to:
  - `color`: Hex code for plot consistency.
  - `label`: Human-readable display name.
  - `direction` (optional): Set to `1` or `-1` to standardize whether a higher score means more or less deleterious.

## 3. config_odds.yaml

Used exclusively by the `top_n_vars_traits` and `top_n_vars_olink` applets. It shares the exact same structure as `config_correlations.yaml` but provides distinct color palettes and label choices optimized for cumulative z-score line plots and heatmaps.
