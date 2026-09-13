"""Source-adapted LDPC components; A003 encoder."""

import math

from a002_qualification import PchkMatrix


def sum_product_decode(matrix, llrs, max_iterations, clip=30.0, early_stop=True):
    """Decode finite LLR evidence with bounded flooding sum-product updates.

    Positive LLR means evidence for bit zero. The returned tuple contains no
    truth or scorer information; parity is only a validity diagnostic.
    """
    if not isinstance(matrix, PchkMatrix):
        raise ValueError("matrix must be a PchkMatrix")
    if (type(matrix.rows) is not int or type(matrix.cols) is not int
            or not 0 < matrix.rows < matrix.cols <= 64800):
        raise ValueError("matrix dimensions must satisfy 0 < rows < cols <= 64800")
    if not isinstance(matrix.entries, tuple):
        raise ValueError("matrix entries must be a tuple")
    if type(max_iterations) is not int or max_iterations < 0:
        raise ValueError("max_iterations must be a nonnegative integer")
    if (type(early_stop) is not bool
            or type(clip) not in (int, float)
            or isinstance(clip, bool)):
        raise ValueError("early_stop must be bool and clip must be numeric")
    try:
        clip = float(clip)
    except (OverflowError, ValueError):
        raise ValueError("clip must be finite and in (0, 30]") from None
    if not math.isfinite(clip) or not 0.0 < clip <= 30.0:
        raise ValueError("clip must be finite and in (0, 30]")
    if not isinstance(llrs, tuple) or len(llrs) != matrix.cols:
        raise ValueError("llrs must be a tuple with one value per variable")
    prior = []
    for value in llrs:
        if isinstance(value, bool) or type(value) not in (int, float):
            raise ValueError("llrs must contain finite numbers")
        try:
            value = float(value)
        except (OverflowError, ValueError):
            raise ValueError("llrs must contain finite numbers") from None
        if not math.isfinite(value):
            raise ValueError("llrs must contain finite numbers")
        prior.append(max(-clip, min(clip, value)))

    check_edges = [[] for _ in range(matrix.rows)]
    variable_edges = [[] for _ in range(matrix.cols)]
    seen = set()
    for edge_index, entry in enumerate(matrix.entries):
        if (not isinstance(entry, tuple) or len(entry) != 2):
            raise ValueError("matrix entries must be pairs")
        row, col = entry
        if (type(row) is not int or type(col) is not int
                or not 0 <= row < matrix.rows or not 0 <= col < matrix.cols
                or entry in seen):
            raise ValueError("matrix entry is invalid or duplicated")
        seen.add(entry)
        check_edges[row].append((edge_index, col))
        variable_edges[col].append(edge_index)

    def clipped(value):
        return max(-clip, min(clip, value))

    def hard_and_syndrome(posterior):
        hard = tuple(0 if value >= 0.0 else 1 for value in posterior)
        syndrome = [0] * matrix.rows
        for row, edges in enumerate(check_edges):
            for _, col in edges:
                syndrome[row] ^= hard[col]
        return hard, tuple(syndrome)

    posterior = tuple(prior)
    hard, syndrome = hard_and_syndrome(posterior)
    syndrome_checks = 1
    if early_stop and not any(syndrome):
        return {"hard_bits": hard, "posterior_llrs": posterior,
                "iterations": 0, "syndrome_weight": 0,
                "parity_satisfied": True, "check_updates": 0,
                "syndrome_checks": syndrome_checks}

    check_to_variable = [0.0] * len(matrix.entries)
    check_updates = 0
    iterations = 0
    for _ in range(max_iterations):
        variable_to_check = [0.0] * len(matrix.entries)
        for col, edges in enumerate(variable_edges):
            total = prior[col]
            for edge_index in edges:
                total += check_to_variable[edge_index]
            for edge_index in edges:
                variable_to_check[edge_index] = clipped(
                    total - check_to_variable[edge_index])

        updated = [0.0] * len(matrix.entries)
        for edges in check_edges:
            degree = len(edges)
            values = [math.tanh(variable_to_check[edge_index] / 2.0)
                      for edge_index, _ in edges]
            prefix = [1.0] * (degree + 1)
            suffix = [1.0] * (degree + 1)
            for index, value in enumerate(values):
                prefix[index + 1] = prefix[index] * value
            for index in range(degree - 1, -1, -1):
                suffix[index] = suffix[index + 1] * values[index]
            for index, (edge_index, _) in enumerate(edges):
                product = prefix[index] * suffix[index + 1]
                if product >= 1.0:
                    message = clip
                elif product <= -1.0:
                    message = -clip
                else:
                    message = clipped(2.0 * math.atanh(product))
                updated[edge_index] = message
        check_to_variable = updated
        posterior = tuple(clipped(
            prior[col] + sum(check_to_variable[edge_index]
                             for edge_index in edges))
            for col, edges in enumerate(variable_edges))
        hard, syndrome = hard_and_syndrome(posterior)
        iterations += 1
        check_updates += len(matrix.entries)
        syndrome_checks += 1
        if early_stop and not any(syndrome):
            break

    weight = sum(syndrome)
    return {"hard_bits": hard, "posterior_llrs": posterior,
            "iterations": iterations, "syndrome_weight": weight,
            "parity_satisfied": weight == 0, "check_updates": check_updates,
            "syndrome_checks": syndrome_checks}


def encode_systematic(matrix: PchkMatrix, payload: tuple[int, ...]) -> tuple[int, ...]:
    """Encode a source suffix using a unit lower-bidiagonal parity prefix."""
    if not isinstance(matrix, PchkMatrix):
        raise ValueError("matrix must be a PchkMatrix")
    if (type(matrix.rows) is not int or type(matrix.cols) is not int
            or not 0 < matrix.rows < matrix.cols <= 64800):
        raise ValueError("matrix dimensions must satisfy 0 < rows < cols <= 64800")
    if (not isinstance(matrix.entries, tuple)
            or any(not isinstance(entry, tuple) or len(entry) != 2
                   for entry in matrix.entries)):
        raise ValueError("matrix entries must be pairs")
    seen = set()
    parity_columns = [[] for _ in range(matrix.rows)]
    payload_columns = [[] for _ in range(matrix.rows)]
    for entry in matrix.entries:
        row, col = entry
        if (type(row) is not int or type(col) is not int
                or not 0 <= row < matrix.rows or not 0 <= col < matrix.cols
                or entry in seen):
            raise ValueError("matrix entry is invalid or duplicated")
        seen.add(entry)
        (parity_columns if col < matrix.rows else payload_columns)[row].append(col)
    for row in range(matrix.rows):
        expected = [row] if row == 0 else [row - 1, row]
        if sorted(parity_columns[row]) != expected:
            raise ValueError("parity block is not unit lower-bidiagonal")
    if (not isinstance(payload, tuple)
            or len(payload) != matrix.cols - matrix.rows
            or any(type(bit) is not int or bit not in (0, 1) for bit in payload)):
        raise ValueError("payload must be a binary tuple of length cols - rows")
    parity = [0] * matrix.rows
    for row in range(matrix.rows):
        value = parity[row - 1] if row else 0
        for col in payload_columns[row]:
            value ^= payload[col - matrix.rows]
        parity[row] = value
    return tuple(parity) + payload
