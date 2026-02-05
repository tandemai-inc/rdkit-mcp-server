# Use an official Python image as the base
ARG BASE_IMAGE=python:3.10-slim
FROM ${BASE_IMAGE}

ARG APP_PORT=8000
ENV DEBIAN_FRONTEND=noninteractive
ARG CN_BUILD=false

WORKDIR /app

# Set up mirrors if CN_BUILD
RUN if [ "$CN_BUILD" = "true" ]; then \
    sed -i 's|http://deb.debian.org/debian|https://mirrors.aliyun.com/debian|g' /etc/apt/sources.list.d/debian.sources && \
    sed -i 's|http://deb.debian.org/debian-security|https://mirrors.aliyun.com/debian-security|g' /etc/apt/sources.list.d/debian.sources; \
    fi
RUN apt-get -y update && apt-get install -y libxrender1 libgl1 libsm6 libxext6 ca-certificates

COPY README.md .
COPY LICENSE .
COPY pyproject.toml .

RUN if [ "$CN_BUILD" = "true" ]; then pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple; fi
RUN pip install .


COPY . .

EXPOSE ${APP_PORT}

# # Set the default command to run your app (update as needed)
CMD ["python", "run_server.py", "--host=0.0.0.0", "--port=8000"]
