# Use an official Python image as the base
FROM public.ecr.aws/docker/library/python:3.10-slim

ARG APP_PORT=8000

WORKDIR /app

# Install git
RUN apt-get -y update && apt-get install -y git libxrender1

COPY README.md .
COPY LICENSE .
COPY pyproject.toml .
RUN pip install .

COPY . .

EXPOSE ${APP_PORT}

# # Set the default command to run your app (update as needed)
CMD ["python", "run_server.py", "--host=0.0.0.0", "--port=8000"]
