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
from primers import (
    reverse_translate,
    design_pcr_primers,
)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024  # 2 MB upload limit


@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    cmp_result = None
    error = None
    sequence = ""
    ref_sequence = ""
    query_sequence = ""
    source_name = ""
    active_tab = "analyze"

    if request.method == "POST":
        action = request.form.get("action", "analyze")
        active_tab = action

        if action == "compare":
            ref_sequence = request.form.get("ref_sequence", "")
            query_sequence = request.form.get("query_sequence", "")
            try:
                if ref_sequence and query_sequence:
                    cmp_result = compare_sequences(ref_sequence, query_sequence)
                else:
                    error = "Please provide both reference and query DNA sequences for comparison."
            except SequenceError as exc:
                error = str(exc)
            except Exception as exc:
                error = f"Comparison failed: {exc}"
        else:
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
        cmp_result=cmp_result,
        error=error,
        sequence=sequence,
        ref_sequence=ref_sequence,
        query_sequence=query_sequence,
        source_name=source_name,
        active_tab=active_tab,
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
        "version": "2.2.0"
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


@app.route("/api/primers", methods=["POST"])
def api_primers():
    """
    REST API endpoint for PCR flanking primer design and Tm analysis.
    Accepts JSON: {"sequence": "ATGG...", "primer_len": 20}
    """
    if not request.is_json:
        return jsonify({"success": False, "error": "Content-Type must be application/json"}), 400

    payload = request.get_json(silent=True) or {}
    raw_seq = payload.get("sequence", "")
    primer_len = int(payload.get("primer_len", 20))

    if not raw_seq:
        return jsonify({"success": False, "error": "Missing 'sequence' field."}), 400

    try:
        report = design_pcr_primers(raw_seq, primer_len=primer_len)
        return jsonify({"success": True, "data": report}), 200
    except ValueError as exc:
        return jsonify({"success": False, "error": str(exc)}), 422
    except Exception as exc:
        return jsonify({"success": False, "error": f"Internal error: {exc}"}), 500


@app.route("/api/reverse-translate", methods=["POST"])
def api_reverse_translate():
    """
    REST API endpoint to reverse-translate amino acid sequence to degenerate DNA.
    Accepts JSON: {"peptide": "MAIV"}
    """
    if not request.is_json:
        return jsonify({"success": False, "error": "Content-Type must be application/json"}), 400

    payload = request.get_json(silent=True) or {}
    peptide = payload.get("peptide", "")

    if not peptide:
        return jsonify({"success": False, "error": "Missing 'peptide' field."}), 400

    try:
        report = reverse_translate(peptide)
        return jsonify({"success": True, "data": report}), 200
    except ValueError as exc:
        return jsonify({"success": False, "error": str(exc)}), 422
    except Exception as exc:
        return jsonify({"success": False, "error": f"Internal error: {exc}"}), 500


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


@app.route("/api/primers", methods=["POST"])
def api_primers():
    """
    REST API endpoint for PCR flanking primer design and Tm analysis.
    Accepts JSON: {"sequence": "ATGG...", "primer_len": 20}
    """
    if not request.is_json:
        return jsonify({"success": False, "error": "Content-Type must be application/json"}), 400

    payload = request.get_json(silent=True) or {}
    raw_seq = payload.get("sequence", "")
    primer_len = int(payload.get("primer_len", 20))

    if not raw_seq:
        return jsonify({"success": False, "error": "Missing 'sequence' field."}), 400

    try:
        report = design_pcr_primers(raw_seq, primer_len=primer_len)
        return jsonify({"success": True, "data": report}), 200
    except ValueError as exc:
        return jsonify({"success": False, "error": str(exc)}), 422
    except Exception as exc:
        return jsonify({"success": False, "error": f"Internal error: {exc}"}), 500


@app.route("/api/reverse-translate", methods=["POST"])
def api_reverse_translate():
    """
    REST API endpoint to reverse-translate amino acid sequence to degenerate DNA.
    Accepts JSON: {"peptide": "MAIV"}
    """
    if not request.is_json:
        return jsonify({"success": False, "error": "Content-Type must be application/json"}), 400

    payload = request.get_json(silent=True) or {}
    peptide = payload.get("peptide", "")

    if not peptide:
        return jsonify({"success": False, "error": "Missing 'peptide' field."}), 400

    try:
        report = reverse_translate(peptide)
        return jsonify({"success": True, "data": report}), 200
    except ValueError as exc:
        return jsonify({"success": False, "error": str(exc)}), 422
    except Exception as exc:
        return jsonify({"success": False, "error": f"Internal error: {exc}"}), 500


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
