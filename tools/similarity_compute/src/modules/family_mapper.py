#
# @Author  : zhuy
# @Date    : 2026/1/10 20:04
# @Desc    : Find all families in a directory
#
from collections import defaultdict
from pathlib import Path
from typing import Dict, List


class FamilyMapper:

    @staticmethod
    def collect_family_dirs(
            root: Path,
            recursive: bool = False,
            return_absolute: bool = False,
            require_underscore: bool = True,
    ) -> Dict[str, List[str]]:
        """
        Scan children under root and group them by family ID.

        Parameters:
          root: Root directory to scan, as a pathlib.Path or a string convertible to Path.
          recursive: Whether to recursively search all subdirectories. Defaults to False.
          return_absolute: Whether values contain absolute paths or only directory names.
          require_underscore: Whether to accept only names containing '_'. If False,
            directories without '_' are treated as entries with key=''.

        Returns:
          dict: { family_id: [folders...] }
        """
        root = Path(root)
        if not root.exists() or not root.is_dir():
            raise ValueError(f"The specified path does not exist or is not a directory: {root}")

        groups = defaultdict(list)

        iterator = root.rglob('*') if recursive else root.iterdir()

        for p in iterator:
            if not p.is_dir():
                continue
            # Use only the subdirectory name. Hidden-directory filters can be added if needed.
            name = p.name.strip()
            if require_underscore:
                if '_' not in name:
                    # Skip directories that do not match <family_id>_<individual>.
                    continue
                family_id, _ = name.split('_', 1)
            else:
                parts = name.split('_', 1)
                family_id = parts[0] if parts else ''
            key = family_id
            val = str(p.resolve()) if return_absolute else name
            groups[key].append(val)

        # Sort by key and sort each value list.
        result = {k: sorted(v) for k, v in sorted(groups.items())}
        return result
