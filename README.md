## PartiotionedPRS custom scripts

### In this repo:
1. Penalized regression models (+ correlation plots)

2. Gene-mapping of identified loci

3. Calculating, goodness-of-fit and association analysis of overall and partitioned PRS

> [!NOTE]
> You need an existing access to UK Biobank genotypic and phenotypic data.
> Preferentially to be run in a Windows env with a WSL (Windows Subsystem for Lunux) installed.
> All dependent Python or R packages should be installed separately

> [!TIP]
> Each phenotype should be a struct with necessary fields: `eid`, `rawUKB`, `numericFlag`, `tag`. Fields `info`, `source`, and `termMeaning` can be left empty. `exeid` field includes exclusion criteria (if any). `eid` for a binary trait contains only ID of cases with that particular trait.
```

UKB_STRUCT_ALL

struct with fields:

           info: [1×1 struct]
            eid: [13884×1 double]
         rawUKB: [13884×1 string]
         source: [13884×1 string]
    numericFlag: 0
          exeid: [906×1 double]
             ex: [1×1 struct]
    termMeaning: [13884×1 string]
            tag: "Btrait"

% for a continuous trait:

UKB_STRUCT_ALL = 

  struct with fields:

           info: [1×1 struct]
            eid: [274350×1 double]
         rawUKB: [274350×1 double]
    numericFlag: 1
            tag: "cTrait"
    termMeaning: ''

```


