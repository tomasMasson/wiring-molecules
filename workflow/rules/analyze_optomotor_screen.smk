import cv2 as cv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from scipy.stats import bootstrap
from pathlib import Path
# Suppress matplotlib warnings
import matplotlib
matplotlib.use('agg')


# Start auxiliar functions


def read_behavioral_data(data_file: str) -> pd.DataFrame:
    column_names = [
                "Frame",
                "Timestamp",
                "Arena",
                "x",
                "y",
                "Orientation",
                "Stimulus",
                "Frequency",
                "Clockwise",
                "Trial",
                "Wedges"]
    return pd.read_csv(data_file, names=column_names, dtype={"x": object, "y": object, "Orientation": object})


def extract_fly_size_from_video(video_file: str):
    # Initialize video in OpenCV
    cap = cv.VideoCapture(video_file)
    # Containers for frame and fly size data
    frame_count = 0
    frames_vector = []
    fly_sizes_vector = []
    # Loop until the end of the video
    while (cap.isOpened()):
        # Capture frame-by-frame
        ret, frame = cap.read()
        # Update frame counter
        frame_count += 1
        # By default, the fly area if NaN
        fly_size = np.nan
        # Finish processing is the video of exhausted
        if frame is None:
            break
        else:
            # Remove upper part of the frame where annotations are
            frame = frame[100:, :]
            # Convert to grayscale
            gray = cv.cvtColor(frame, cv.COLOR_BGR2GRAY)
            # Threshold image
            _, dst = cv.threshold(gray, 110, 255, cv.THRESH_BINARY)
            # Select contour compatible with fly size
            contours, _ = cv.findContours(dst, cv.RETR_TREE, cv.CHAIN_APPROX_NONE)
            for contour in contours:
                if len(contour) > 10 and len(contour) < 200:
                    # Get fly contour area
                    fly_size = cv.contourArea(contour)
        frames_vector.append(frame_count)
        fly_sizes_vector.append(fly_size)
    # release the video capture object
    cap.release()
    # Closes all the windows currently opened.
    cv.destroyAllWindows()
    # Save results as output file
    df = pd.DataFrame({"Video_frame": frames_vector, "Size": fly_sizes_vector})
    df["Fly_number"] = Path(video_file).name.split("_")[0]

    return df


def add_fly_size_to_table(df: pd.DataFrame, data_file: str) -> pd.DataFrame:
    video_file = "resources/optomotor_screen_dataset/videos/" + Path(data_file).name[:-3] + "mp4"
    sdf = extract_fly_size_from_video(video_file)
    return (df
            .assign(Fly_size=sdf["Size"],
                    Video_frame=sdf["Video_frame"]
                    )
            )


def assign_fly_number(df: pd.DataFrame, data_file: str) -> pd.DataFrame:
    return (df
            .assign(Fly_number=Path(data_file).name.split("_")[0])
            )


def assign_trial_id(df: pd.DataFrame) -> pd.DataFrame:
    return (df
            .assign(Trial_id=df["Fly_number"].astype("str") + "_" + df["Trial"].astype("str"))
            )


def assign_perturbation_id(df: pd.DataFrame, perturbations_file: str) -> pd.DataFrame:
    pdf = pd.read_csv(perturbations_file)
    return (df
            .assign(Perturbation=df.astype({"Fly_number": "int32"}).Fly_number.replace(dict(zip(pdf["Fly_number"], pdf["Genotype"]))))
            )


def interpolate_missing_values(df: pd.DataFrame) -> pd.DataFrame:
    # Missing values are labeled as "None" in raw data
    # Maximum window of Nan values is ten frames
    return df.replace({"None": np.nan}).astype({"x": "float64", "y": "float64", "Orientation": "float64"}).interpolate(limit=10, axis=0)


def drop_missing_data_trials_1(df: pd.DataFrame) -> pd.DataFrame:
    # Remove trials where 10 frames or more are missing during stimulation
    trials_to_drop = [trial
                      for trial in df["Trial"].unique()
                      if sum(df[(df["Trial"] == trial)&(df["Stimulus"] == "on")].x.isna()) >= 10]
    return df[~df["Trial"].isin(trials_to_drop)]


def drop_missing_data_trials_2(df: pd.DataFrame) -> pd.DataFrame:
    # Remove trials where 60 frames or more are missing during pre-stimulation
    trials_to_drop = [trial
                      for trial in df["Trial"].unique()
                      if sum(df[(df["Trial"] == trial)&(df["Stimulus"] == "off")].x.isna()) >= 60]
    return df[~df["Trial"].isin(trials_to_drop)]


def drop_close_to_wall_trials(df: pd.DataFrame) -> pd.DataFrame:
    df["Distance_center"] = np.sqrt(df["x"]**2 + df["y"]**2)
    trials_to_drop = [trial
                      for trial in df["Trial"].unique()
                      if sum(df[(df["Trial"] == trial)]["Distance_center"] > 25) > 30]
    return df[~df["Trial"].isin(trials_to_drop)]
    return df


def transform_pixel_to_mm(df: pd.DataFrame, pixel_resolution: int) -> pd.DataFrame:
    return (df
            .assign(x=df.x / pixel_resolution,
                    y=df.y / pixel_resolution
                   )
           )


def correct_turning_response(df: pd.DataFrame) -> pd.DataFrame:
    # Rectifies head-rear ambiguity from OpenCV
    turning_responses = []
    for difference in df.Orientation.diff():
        if difference > 90:
            turning_responses.append(difference - 180)
        elif difference < -90:
            turning_responses.append(180 + difference)
        else:
            turning_responses.append(difference)
    return (df
            .assign(Turning_responses=turning_responses)
           )


def correct_responses_by_pinwheel_sense(df: pd.DataFrame) -> pd.DataFrame:
    # Invert x coordinate and turning responses when pinwheel is counter-clockwise
    return (df
            .assign(x=df.x * df.Clockwise.replace({False: -1, True: 1}),
                    Turning_responses=df.Turning_responses * df.Clockwise.replace({False: -1, True: 1})
                    )
            )


def compute_partial_speeds(df: pd.DataFrame) -> pd.DataFrame:
    return (df
            .assign(dt=df.Timestamp.diff() * 1e9,
                    vx=(df.x.diff() / df.Timestamp.diff() * 1e9),
                    vy=(df.y.diff() / df.Timestamp.diff() * 1e9)
                   )
           )


def compute_ground_speed(df: pd.DataFrame) -> pd.DataFrame:
    return (df
            .assign(Ground_speed=np.sqrt(df.vx**2 + df.vy**2))
           )


def drop_inactive_trials(df: pd.DataFrame) -> pd.DataFrame:
    # Keep trials where the average ground speed is at least 5 mm per second
    return df.groupby(["Trial_id"]).filter(lambda x: x["Ground_speed"].mean() >= 5)


def drop_inactive_flies(df: pd.DataFrame) -> pd.DataFrame:
    # Keep experimental flies where at least 10 trials are active
    active_flies = list(df["Fly_number"].value_counts().loc[lambda x: x >= 3600].index)
    return df[df["Fly_number"].isin(active_flies)]


def compute_cumulative_turning(df: pd.DataFrame) -> pd.DataFrame:
    dff = pd.DataFrame()
    for trial in df["Trial"].unique():
        tdf = df[df["Trial"] == trial]
        # Trim all trials to 360 frames (2s pre + 4s stimulus)
        tdf = tdf.iloc[-360:, :]
        tdf["Frame"] = list(range(1, 361))
        # Compute cumulative turning response (aligned to stimulus onset)
        tdf["Cumulative_turning"] = tdf["Turning_responses"].cumsum() - tdf["Turning_responses"].cumsum().iloc[120]
        dff = pd.concat([dff, tdf])
    return (dff)


def process_data_file(data_file: str) -> pd.DataFrame:
  return (read_behavioral_data(data_file)
          .pipe(add_fly_size_to_table, data_file)
          .pipe(assign_fly_number, data_file)
          .pipe(assign_trial_id)
          .pipe(assign_perturbation_id, "resources/fly_experimental_labels.csv")
          .pipe(interpolate_missing_values)
          .pipe(drop_missing_data_trials_1)
          .pipe(drop_missing_data_trials_2)
          .pipe(transform_pixel_to_mm, 9.5)
          .pipe(drop_close_to_wall_trials)
          .pipe(correct_turning_response)
          .pipe(correct_responses_by_pinwheel_sense)
          .pipe(compute_partial_speeds)
          .pipe(compute_ground_speed)
          .pipe(drop_inactive_trials)
          .pipe(drop_inactive_flies)
          .pipe(compute_cumulative_turning)
          )


def pool_folder_data(folder: str) -> pd.DataFrame:
    pooled_df = pd.DataFrame()
    for file in Path(folder).iterdir():
        pooled_df = pd.concat([pooled_df, process_data_file(file)])
    return pooled_df


def save_processed_data(folder: str, output_file:str):
    pool_folder_data(folder).to_csv(output_file, index=False)


def pool_dataframes(file_list: list) -> pd.DataFrame:
    df_list = [pd.read_csv(file) for file in file_list]
    return pd.concat(df_list)


def compute_optomotor_index(array: pd.Series) -> pd.Series:
    return array.sum() / sum(array.abs())


def plot_optomotor_index(file_list: list, output_plot: str):
    df = pool_dataframes(file_list)
    sns.set(style="ticks", rc={"font.family": "Carlito", 'axes.linewidth': 1})
    fig, axs = plt.subplots(figsize=(4,2))
    # Optomotor index
    trial_avg = df[df["Stimulus"] == "on"].groupby(["Perturbation", "Fly_number", "Trial"])["Turning_responses"].apply(compute_optomotor_index).reset_index(name="Optoindex")
    fly_avg = trial_avg.groupby(["Perturbation", "Fly_number"])["Optoindex"].mean().reset_index()
    plot_order = list(fly_avg.groupby(["Perturbation"])["Optoindex"].mean().sort_values(ascending=False).index)
    print(plot_order)
    sns.stripplot(x="Perturbation",
                  y="Optoindex",
                  data=fly_avg,
                  order=plot_order,
                  color="black",
                  alpha=.4,
                  ax=axs)
    sns.pointplot(x="Perturbation",
                  y="Optoindex",
                  data=fly_avg,
                  order=plot_order,
                  color="black",
                  linestyle="none",
                  ax=axs)
    sns.despine(bottom=True)
    plt.xlabel("")
    plt.tick_params(axis="x", which="major", width=0, labelsize=0, rotation=60)
    plt.ylim(-0.05, 1)
    plt.ylabel("Optomotor index", fontsize=8)
    plt.tick_params(axis="y", which="major", width=1, labelsize=8)
    plt.tight_layout()
    fig.savefig(output_plot)


def plot_angular_speed(file_list: list, output_plot: str):
    df = pool_dataframes(file_list)
    sns.set(style="ticks", rc={"font.family": "Carlito", 'axes.linewidth': 1})
    fig, axs = plt.subplots(figsize=(4,2))
    # Optomotor index
    trial_avg = df[df["Stimulus"] == "on"].groupby(["Perturbation", "Fly_number", "Trial"])["Cumulative_turning"].max().reset_index(name="Angular_speed")
    fly_avg = trial_avg.groupby(["Perturbation", "Fly_number"])["Angular_speed"].mean().reset_index()
    fly_avg["Angular_speed"] = fly_avg["Angular_speed"] / 4
    plot_order = list(fly_avg.groupby(["Perturbation"])["Angular_speed"].mean().sort_values(ascending=False).index)
    sns.stripplot(x="Perturbation",
                  y="Angular_speed",
                  data=fly_avg,
                  order=plot_order,
                  color="black",
                  alpha=.4,
                  ax=axs)
    sns.pointplot(x="Perturbation",
                  y="Angular_speed",
                  data=fly_avg,
                  order=plot_order,
                  color="black",
                  linestyle="none",
                  ax=axs)
    sns.despine(bottom=True)
    plt.xlabel("")
    plt.xticks(ticks=[], labels=[])
    plt.tick_params(axis="x", which="major", width=0, labelsize=0, rotation=60)
    plt.ylim(0, 400)
    plt.ylabel("Angular speed", fontsize=8)
    plt.tick_params(axis="y", which="major", width=1, labelsize=8)
    plt.tight_layout()
    fig.savefig(output_plot)


def plot_response_threshold_frame(input_file: str, output_plot: str):
    df = pd.read_csv(input_file)
    sns.set(style="ticks", rc={"font.family": "Carlito", "figure.figsize": (3, 2), 'axes.linewidth': 1})
    fig, axs = plt.subplots(ncols=1)
    # Optomotor index
    trial_avg = df[(df["Frame"] > 120)&(df["Cumulative_turning"] > 200)].groupby(["Perturbation", "Fly_number", "Trial"])["Frame"].min().reset_index(name="Threshold")
    fly_avg = trial_avg.groupby(["Perturbation", "Fly_number"])["Threshold"].mean().reset_index()
    fly_avg["Threshold"] = fly_avg["Threshold"] * (1/60) - 2
    sns.stripplot(x="Perturbation",
                  y="Threshold",
                  data=fly_avg,
                  color="black",
                  alpha=.4,
                  ax=axs)
    sns.pointplot(x="Perturbation",
                  y="Threshold",
                  data=fly_avg,
                  color="black",
                  ax=axs)
    sns.despine(ax=axs)
    axs.set_ylim(0, 4)
    plt.tight_layout()
    fig.savefig(output_plot)


def plot_anti_optomotor_response(input_file: str, output_plot: str):
    df = pd.read_csv(input_file)
    sns.set(style="ticks", rc={"font.family": "Carlito", "figure.figsize": (3, 2), 'axes.linewidth': 1})
    fig, axs = plt.subplots(ncols=1)
    # Optomotor index
    trial_avg = df[(df["Stimulus"] == "on")&(df["Turning_responses"] < 0)].groupby(["Perturbation", "Fly_number", "Trial"])["Turning_responses"].sum().reset_index()
    fly_avg = trial_avg.groupby(["Perturbation", "Fly_number"])["Turning_responses"].mean().reset_index()
    sns.stripplot(x="Perturbation",
                  y="Turning_responses",
                  data=fly_avg,
                  color="black",
                  alpha=.4,
                  ax=axs)
    sns.pointplot(x="Perturbation",
                  y="Turning_responses",
                  data=fly_avg,
                  color="black",
                  ax=axs)
    sns.despine(ax=axs)
    axs.set_ylim(0, -1000)
    plt.tight_layout()
    fig.savefig(output_plot)


def extract_saccades(data_table: str, output_file: str):
    df = pd.read_csv(data_table)
    # Saccades table
    sdf = pd.DataFrame()
    for trial_id in df["Trial_id"].unique():
        tdf = df[(df["Trial_id"] == trial_id)]
        # Get syn saccades
        syn_peaks = find_peaks(savgol_filter(tdf["Turning_responses"], 21, 4), prominence=1, height=5, width=3)
        # Get anti saccades
        anti_peaks = find_peaks(-1 * savgol_filter(tdf["Turning_responses"], 21, 4), prominence=1, height=5, width=3)
        sdf = pd.concat([sdf, pd.DataFrame({"Center": syn_peaks[0], "Height": syn_peaks[1]["peak_heights"], "Type": ["Syn"] * len(syn_peaks[0]), "Trial_id": [trial_id] * len(syn_peaks[0])}), pd.DataFrame({"Center": anti_peaks[0], "Height": anti_peaks[1]["peak_heights"], "Type": ["Anti"] * len(anti_peaks[0]), "Trial_id": [trial_id] * len(anti_peaks[0])})])
    sdf.to_csv(output_file, index=False)


# End auxiliar functions


# PERTURBATIONS = ["Urbana_wt", "mCherry_control", "Fas2_RNAi", "Kek1_RNAi", "Lac_RNAi", "Lac_UAS", "Ptp10D_RNAi", "Ptp99A_RNAi_39006", "Ptp99A_RNAi_57299", "Sdc_RNAi", "Con_UAS", "Toll-7_UAS", "Kek1_Ptp10D_RNAix2", "Kek1_RNAi+Con_UAS", "Kek1_RNAi+Elff_UAS", "Kek1_Nlg3_RNAix2+Con_UAS", "Kek1_Fas2_Ptp99a_RNAix3", "T4T5_Kir2.1_R42F06"]
PERTURBATIONS = ["Urbana_wt", "mCherry_control", "Kek1_Ptp10D_RNAix2"]


rule all:
    input:
        "results/optomotor_behavioral_screen/screen_optomotor_index.svg",
        "results/optomotor_behavioral_screen/screen_turning_response.svg",
        [f"results/optomotor_behavioral_screen/{perturbation}_anti_optomotor_response.svg" for perturbation in PERTURBATIONS],
        [f"results/optomotor_behavioral_screen/{perturbation}_saccades.csv" for perturbation in PERTURBATIONS]


rule preprocess_behavioral_data:
    input:
        [f"resources/optomotor_screen_dataset/{perturbation}" for perturbation in PERTURBATIONS]
    output:
        [f"results/optomotor_behavioral_screen/{perturbation}_pooled.csv" for perturbation in PERTURBATIONS]
    run:
        for inp, out in zip(input, output):
            save_processed_data(inp, out)


rule plot_optomotor_index:
    input:
        [f"results/optomotor_behavioral_screen/{perturbation}_pooled.csv" for perturbation in PERTURBATIONS]
    output:
        "results/optomotor_behavioral_screen/screen_optomotor_index.svg"
    run:
        plot_optomotor_index(input, output[0])


rule plot_angular_speed:
    input:
        [f"results/optomotor_behavioral_screen/{perturbation}_pooled.csv" for perturbation in PERTURBATIONS]
    output:
        "results/optomotor_behavioral_screen/screen_turning_response.svg"
    run:
        plot_angular_speed(input, output[0])


rule plot_response_threshold:
    input:
        [f"results/optomotor_behavioral_screen/{perturbation}_pooled.csv" for perturbation in PERTURBATIONS]
    output:
        [f"results/optomotor_behavioral_screen/{perturbation}_response_threshold.svg" for perturbation in PERTURBATIONS]
    run:
        for inp, out in zip(input, output):
            plot_response_threshold_frame(inp, out)


rule plot_anti_optomotor_response:
    input:
        [f"results/optomotor_behavioral_screen/{perturbation}_pooled.csv" for perturbation in PERTURBATIONS]
    output:
        [f"results/optomotor_behavioral_screen/{perturbation}_anti_optomotor_response.svg" for perturbation in PERTURBATIONS]
    run:
        for inp, out in zip(input, output):
            plot_anti_optomotor_response(inp, out)


rule extract_saccades:
    input:
        [f"{perturbation}_pooled.csv" for perturbation in PERTURBATIONS]
    output:
        [f"results/optomotor_behavioral_screen/{perturbation}_saccades.csv" for perturbation in PERTURBATIONS]
    run:
        for inp, out in zip(input, output):
            extract_saccades(inp, out)
