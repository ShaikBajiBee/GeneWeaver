import json
from io import BytesIO
from flask import (
    Flask,
    render_template,
    request,
    jsonify,
    Response,
    send_file,
)
from analyzer import (
    analyze_sequence,
    parse_fasta,
    generate_text_report,
    compare_sequences,
    SequenceError,
)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024  # 2 MB upload limit


@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    error = None
    sequence = ""
    source_name = ""
    records = []
    selected_record_idx = 0

    if request.method == "POST":
        sequence = request.form.get("sequence", "")
        uploaded = request.files.get("fasta_file")

        try:
            if uploaded and uploaded.filename:
                source_name = uploaded.filename
                content = uploaded.read().decode("utf-8-sig")
                records = parse_fasta(content)
                if records:
                    sequence = records[0]["sequence"]
                    source_name = f"{uploaded.filename} ({records[0]['header']})"
            
            if sequence:
                result = analyze_sequence(sequence)
            else:
                error = "Please provide a DNA sequence or upload a file."
        except UnicodeDecodeError:
            error = "The uploaded file must be a valid UTF-8 text or FASTA file."
        except SequenceError as exc:
            error = str(exc)
        except Exception as exc:
            error = f"Could not process input: {exc}"

    return render_template(
        "index.html",
        result=result,
        error=error,
        sequence=sequence,
        source_name=source_name,
    )


# -------------------------------------------------------------
# REST API Endpoints
# -------------------------------------------------------------
@app.route("/api/health", methods=["GET"])
def api_health():
    """Health check endpoint for monitoring."""
    return jsonify({
        "status": "healthy",
        "service": "GeneWeaver DNA Analysis API",
        "version": "2.0.0"
    }), 200


@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    """
    REST API endpoint for programmatic DNA sequence analysis.
    Accepts JSON: {"sequence": "ATGGCC..."}
    """
    if not request.is_json:
        return jsonify({
            "success": False,
            "error": "Request must have Content-Type: application/json"
        }), 400

    payload = request.get_json(silent=True) or {}
    raw_sequence = payload.get("sequence", "")

    if not raw_sequence:
        return jsonify({
            "success": False,
            "error": "Missing 'sequence' field in JSON payload."
        }), 400

    try:
        report = analyze_sequence(raw_sequence)
        return jsonify({
            "success": True,
            "data": report
        }), 200
    except SequenceError as exc:
        return jsonify({
            "success": False,
            "error": str(exc)
        }), 422
    except Exception as exc:
        return jsonify({
            "success": False,
            "error": f"Internal server error: {exc}"
        }), 500


@app.route("/api/compare", methods=["POST"])
def api_compare():
    """
    REST API endpoint for sequence alignment & variant/mutation detection.
    Accepts JSON: {"reference": "ATGG...", "query": "ATGC..."}
    """
    if not request.is_json:
        return jsonify({
            "success": False,
            "error": "Request must have Content-Type: application/json"
        }), 400

    payload = request.get_json(silent=True) or {}
    ref = payload.get("reference", "")
    query = payload.get("query", "")

    if not ref or not query:
        return jsonify({
            "success": False,
            "error": "Both 'reference' and 'query' sequence fields are required."
        }), 400

    try:
        report = compare_sequences(ref, query)
        return jsonify({
            "success": True,
            "data": report
        }), 200
    except SequenceError as exc:
        return jsonify({
            "success": False,
            "error": str(exc)
        }), 422
    except Exception as exc:
        return jsonify({
            "success": False,
            "error": f"Internal server error: {exc}"
        }), 500


# -------------------------------------------------------------
# Export Endpoints
# -------------------------------------------------------------
@app.route("/export/json", methods=["POST"])
def export_json():
    """Download analysis report as JSON."""
    sequence = request.form.get("sequence", "")
    try:
        report = analyze_sequence(sequence)
        buffer = BytesIO(json.dumps(report, indent=2).encode("utf-8"))
        return send_file(
            buffer,
            mimetype="application/json",
            as_attachment=True,
            download_name="geneweaver_report.json",
        )
    except Exception as exc:
        return Response(f"Export failed: {exc}", status=400)


@app.route("/export/txt", methods=["POST"])
def export_txt():
    """Download analysis report as formatted text."""
    sequence = request.form.get("sequence", "")
    try:
        report = analyze_sequence(sequence)
        text_content = generate_text_report(report)
        buffer = BytesIO(text_content.encode("utf-8"))
        return send_file(
            buffer,
            mimetype="text/plain",
            as_attachment=True,
            download_name="geneweaver_summary.txt",
        )
    except Exception as exc:
        return Response(f"Export failed: {exc}", status=400)


@app.errorhandler(413)
def too_large(_error):
    return jsonify({"error": "File exceeds maximum 2 MB limit."}), 413


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
