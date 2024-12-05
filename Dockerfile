FROM python:3.10

# Base directory
ENV docker_work_dir=/app
WORKDIR /app

# Data directory
ENV docker_data_dir=/data

# Input and output directories
ENV docker_input_dir=/data/updated_mmcifs

# Make input directory
RUN mkdir -p /data/output/updated_mmcifs

# Copy files
COPY ./requirements.txt ${docker_work_dir}
COPY ./find_conformers.py ${docker_work_dir}
COPY cluster_conformers ${docker_work_dir}/cluster_conformers

RUN pip install wheel
RUN pip install --no-cache-dir -r requirements.txt --verbose

ENTRYPOINT [ "python", "find_conformers.py", "-c",  "/data/output/ca_distances", "-d", "/data/output/distance_differences", "-s", "/data/output/cluster_results"  ]
CMD [ ]

# Example:
#CMD [ "-u", "O34926", "-m", "/data/updated_mmcifs/3nc3_updated.cif", "A", "B", "-m", "/data/updated_mmcifs/3nc5_updated.cif", "A", "B", "-m", "/data/updated_mmcifs/3nc6_updated.cif", "A", "B", "-m", "/data/updated_mmcifs/3nc7_updated.cif", "A", "B"  ]