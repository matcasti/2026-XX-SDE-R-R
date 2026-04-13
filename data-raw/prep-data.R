
library(data.table)
library(RHRV)
library(CardioCurveR)

subjects <- list.files("data-raw/fantasia/") |>
  gsub(pattern = "\\.hea|\\.dat|\\.ecg|-", replacement = "") |>
  unique()

for (id in subjects) {
  rr_intervals <- RHRV::LoadBeatWFDB(
    HRVData = RHRV::CreateHRVData(),
    RecordName = paste0("data-raw/fantasia/",id),
    annotator = "ecg"
  )$Beat$Time |>
    diff() |> CardioCurveR::clean_outlier(threshold = 5)

cat(rr_intervals, sep = "\n", file = paste0("data/fantasia/",id,".txt"))
}
