class GlobUtils {
    // Convert glob pattern to Java regex.
    // Escapes all regex meta-characters first, then restores * and ? as wildcards.
    static String globToRegex(String glob) {
        def escaped = glob
            .replace('\\', '\\\\')  // must be first
            .replace('.', '\\.')
            .replace('+', '\\+')
            .replace('(', '\\(')
            .replace(')', '\\)')
            .replace('[', '\\[')
            .replace(']', '\\]')
            .replace('{', '\\{')
            .replace('}', '\\}')
            .replace('^', '\\^')
            .replace('$', '\\$')
            .replace('|', '\\|')
            .replace('?', '.')
            .replace('*', '.*')
        return /^${escaped}$/
    }
}
