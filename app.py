from flask import Flask, render_template, request
from analyzer import analyze_sequence, parse_fasta, SequenceError

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024  # 2 MB upload limit

@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    error = None
    sequence = ""
    source_name = ""

    if request.method == "POST":
        sequence = request.form.get("sequence", "")
        uploaded = request.files.get("fasta_file")

        try:
            if uploaded and uploaded.filename:
                source_name = uploaded.filename
                content = uploaded.read().decode("utf-8-sig")
                sequence = parse_fasta(content)
            result = analyze_sequence(sequence)
        except UnicodeDecodeError:
            error = "The uploaded file must be a UTF-8 text FASTA file."
        except SequenceError as exc:
            error = str(exc)
        except Exception as exc:
            error = f"Could not process the input: {exc}"

    return render_template(
        "index.html",
        result=result,
        error=error,
        sequence=sequence,
        source_name=source_name,
    )

@app.errorhandler(413)
def too_large(_error):
    return "File too large. Please upload a file smaller than 2 MB.", 413

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
